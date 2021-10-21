#! /bin/tcsh -f
#
#  add masked version of Fo-Fc map to partial structure of the solvent
#
#
set pdbfile = ""
set mapfile = ""
set mtzfile = ""
set lastlog = ""
set outprefix = "touchedup"

set diffmapfile = fofc_small.map

set trials = 20
set scaleoverall = 1
set Boverall = 0
set min_delta_rho = 0.0001
set keepmaps = 0

# smooth the fofc map
set fftB = 0

set mask_scale = 1
set maskB = 0
set mask_reso = 1.3

set Bramp = 0
set iterate = 200
set skiprefine = 0

# define native Shannon voxel size (median atomic B factor)
set solventB = 5

set super_mult = ""

set tempfile = /dev/shm/${USER}/touchupsolvent_$$_
mkdir -p /dev/shm/${USER}

set do_ramp = 0
set damp_ramp = 0

set logfile = touchupsolvent_details.log

echo "command-line: $* "

foreach arg ( $* )
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set Val = `echo $arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if("$Key" =~ *.pdb && ! $assign )  set pdbfile = $Key
    if("$key" == "pdbfile") set pdbfile = "$Val"

    if("$Key" =~ *.mtz && ! $assign)  set mtzfile = $Key
    if("$key" == "mtzfile") set mtzfile = "$Val"
    if("$Key" =~ *.map && ! $assign)  set mapfile = $Key
    if("$key" == "mapfile") set mapfile = "$Val"

    if("$key" =~ outpre* && $assign) set outprefix = "$Val"

    if("$key" == "trials" && $assign) set trials = "$num"
    if("$key" == "scale" && $assign) set scaleoverall = "$num"
    if("$Key" == "B" && $assign) set Boverall = "$num"
    if("$Key" == "Boverall" && $assign) set Boverall = "$num"

    if("$Key" == "fftB" && $assign) set fftB = "$num"
    if("$Key" == "fft_B" && $assign) set fftB = "$num"

    if("$Key" == "mask_scale" && $assign) set mask_scale = "$num"
    if("$Key" == "maskB" && $assign) set maskB = "$num"
    if("$Key" == "mask_B" && $assign) set maskB = "$num"
    if("$Key" == "mask_reso" && $assign) set mask_reso = "$num"

    if("$key" == "bramp" && "$int" != "0" ) set Bramp = "$int"

    if("$key" == "iterate" && "$int" != "0" ) set iterate = "$int"
    if("$key" == "iterate" && "$val" == "ramp") set do_ramp = 1
    if("$key" == "iterate" && "$val" == "damp") set damp_ramp = 1
    if("$key" == "iterate" && "$val" == "skip") set skiprefine = 1
    if("$key" == "skiprefine") set skiprefine = 1
    if("$key" == "nramp") set do_ramp = "$num"
    if("$key" == "damp") set damp_ramp = "$num"
    if("$key" == "keepmaps") set keepmaps = "$num"


    if("$key" == "tempfile" && $assign) then
        set tempfile = "$Val"
        setenv DEBUG
    endif
    if("$key" == "debug") then
        setenv DEBUG
    endif
end

if(! -e "$pdbfile" && "$pdbfile" != "") then
    echo "WARNING: $pdbfile does not exist"
endif
if(! -e "$pdbfile") then
    set pdbfile = `ls -1rt *.pdb | egrep "refmacout|minimized" | grep -v _small.pdb | tail -n 1`
    echo "selected $pdbfile"
endif
if(! -e "$mapfile") then
    set mapfile = `ls -1rt *.map | egrep "^total" | tail -n 1`
    echo "selected $mapfile"
endif
if(! -e "$pdbfile") then
    set BAD = "no pdb file"
    goto exit
endif
if(! -e "$mapfile") then
    set BAD = "no map file"
    goto exit
endif


set Bsched
if( "$Bramp" != "0" ) then
    # user-specified B-factor list
    set Bsched = `echo $Bramp | awk '{gsub(","," ");print}'`
endif
if( $#Bsched <= 1 ) then
    # generate a list automatically, if user just said, "1", make it a 50-step ramp
    if( $Bramp <= 1 ) set Bramp = 50
    # generate a list of n=Bramp B factors
    set Bsched = `echo $Bramp | awk '$1>1{step=exp(log(10000./2)/($1-1));for(B=10000;B>2-1e-9;B/=step)print sprintf("%.2g\n",B)+0}'`
    echo "default B factor schedule: $Bsched"
endif
if( $#Bsched >= 1 ) then
    set iterate = $#Bsched
endif


if(! -e refmac_opts.txt) then
  cat << EOF >! refmac_opts.txt
solvent no
SCPART 1
weight matrix 1
bfactor set 5
bfactor 10
damp 0.1 1e-6
#monitor vdw 3
#extern harmonic residues from 1 A to 3072 A sigma 0.02
EOF
  echo "created refmac_opts.txt"
endif


# figure out supercell
set CELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`

echo | mapdump mapin $mapfile >! ${tempfile}mapdump.txt

set GRID = `awk '/Grid sampling/{print "GRID", $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump.txt`
set smallCELL = `awk '/Cell dimensions /{print $4,$5,$6,$7,$8,$9;exit}' ${tempfile}mapdump.txt`
set smallSGnum  = `awk '/ Space-group /{print $NF;exit}' ${tempfile}mapdump.txt`
set smallSG = `awk -v n=$smallSGnum '$1==n{print $4;exit}' ${CLIBD}/symop.lib`

if(! $?super_mult) set super_mult
if(! $?smallCELL && $#super_mult == 3 ) then
    set smallCELL = `echo $CELL $super_mult | awk '{print $1/$7,$2/$8,$3/$9,$4,$5,$6}'`
endif
if(! $?smallCELL ) then
    set smallCELL = ( $CELL )
endif

if(! $?CELL && $#smallCELL == 6 && $#super_mult == 3) then
    set CELL = `echo $smallCELL $super_mult | awk '{print $1*$7,$2*$8,$3*$9,$4,$5,$6}'`
endif
if( $#super_mult != 3 ) then
    set super_mult = `echo $CELL $smallCELL | awk '{printf("%.0f %.0f %.0f", $1/$7,$2/$8,$3/$9)}'`
endif
if(! $?smallSGnum) then
    set smallSGnum = `awk -v sg=$smallSG '$4==sg{print $1;exit}' ${CLIBD}/symop.lib`
endif
set smallSG = `awk -v n=$smallSGnum '$1==n{print $4;exit}' ${CLIBD}/symop.lib`
set ncells = `echo $super_mult | awk '{printf("%.0f\n",$1*$2*$3)}'`

echo $smallCELL |\
awk 'NF==6{DTR=atan2(1,1)/45; A=cos(DTR*$4); B=cos(DTR*$5); G=cos(DTR*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf "%.3f\n", $1*$2*$3*sqrt(skew)}' |\
cat >! ${tempfile}volume
set smallCELLvolume = `cat ${tempfile}volume`
rm -f ${tempfile}volume
set CELLvolume = `echo $smallCELLvolume $ncells | awk '{print $1*$2}'`

set symops = `awk -v n=$smallSGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`
set nasus = `echo $ncells $symops | awk '{print $1*$2}'`

set pdbSG = `awk -v SGnum="$smallSGnum" -F "[\047]" '$1+0==SGnum{print $2;exit}' ${CLIBD}/symop.lib`
echo $smallCELL $pdbSG |\
awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' |\
cat >! smallcell.pdb
echo $CELL |\
awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n",$1,$2,$3,$4,$5,$6)}' |\
cat >! bigcell.pdb



set n = `ls -1rt start[0-9]*.pdb | awk '{print substr($1,match($1,"[0-9]"))+0}' | sort -g | tail -n 1`
if( "$n" == "") set n = 0
#set s = `ls -1rt total_solvent[0-9]*.map | awk '{print substr($1,match($1,"[0-9]"))+0}' | sort -g | tail -n 1`
if( $skiprefine ) goto converged

echo "starting with n= $n"

echo "$pdbfile -> start${n}.pdb"
cp $pdbfile ${tempfile}.pdb
cp ${tempfile}.pdb start${n}.pdb
echo "$mapfile -> total_solvent${n}.map"
cp $mapfile ${tempfile}.map
cp ${tempfile}.map total_solvent${n}.map


# refine first to make sure map is on right scale
set cycle = 0
while ( $cycle <= $iterate )

@ cycle = ( $cycle + 1 )

if( "$Bramp" != "0" ) then
 if( $#Bsched >= $cycle ) then
   set B = $Bsched[$cycle]
   set fftB = $B
   set maskB = $B
   echo "fftB = $fftB"
   echo "maskB = $maskB"
 else
   echo "B factor ramp ended"
   set Bramp = 0
 endif
endif

set n = `ls -1rt start*.pdb | tail -n 1 | awk '{print substr($0,6)+0}'`
echo "n= $n"

set scaled = 0
foreach itr ( `seq 1 10` )

  # transform from master solvent file
  map_scaleB_runme.com $mapfile scale=$scaleoverall B=$Boverall output=total_solvent.map
  solmap2refme.com total_solvent.map start${n}.pdb

  if( $damp_ramp ) then
    damp_ramp_runme.com trials=$trials | tee damp_ramp${n}.log
    set n = `ls -1 start*.pdb converge*.log | awk '{print substr($1,match($1,"[0-9]"))+0}' | sort -g | tail -n 1`
    cp refmacout.pdb start${n}.pdb
  endif

  if( $do_ramp ) then
    ramp_runme.com start${n}.pdb trials=10 occreject=0 \
      weight=1 maxweight=10 minweight=0.1 weightstep=2 |\
    tee ramp${n}.log
    set n = `ls -1 start*.pdb converge*.log | awk '{print substr($1,match($1,"[0-9]"))+0}' | sort -g | tail -n 1`
    cp refmacout.pdb start${n}.pdb
  endif

  set n = `ls -1rt start*.pdb | tail -n 1 | awk '{print substr($0,6)+0}'`
  echo "converge $n"
  converge_refmac.com start${n}.pdb refme.mtz append nosalvage trials=$trials >! converge${n}.log
  tail -q -n 1 refmac_Rplot.txt refmac_shifts.txt refmac_scales.log

  set scale = `tail -n 1 refmac_scales.log | awk 'NF>=4{print $4+0}'`
  set B = `tail -n 1 refmac_scales.log | awk 'NF>=5{print $5+0}'`
  if( "$scale" == "" ) set scale = 1
  if( "$B" == "" ) set B = 0
  set scaleoverall = `echo $scaleoverall $scale | awk '{print $1*$2}'`
  set Boverall = `echo $Boverall $B | awk '{print $1+$2}'`
  echo "$n $scaleoverall $Boverall" | tee -a overall_scaleB.txt

  @ n = ( $n + 1 )
  cp refmacout.pdb start${n}.pdb
  set reindex = `echo $super_mult | awk '{print "h/"$1",k/"$2",l/"$3}'`
  reindex hklin refmacout.mtz hklout start${n}_small.mtz << EOF >! $logfile
reindex $reindex
symm $smallSGnum
EOF


  set test = `echo $scale $B | awk '{print ( sqrt(($1-1)^2)<0.01 && sqrt(($2)^2)<0.1 ) }'`
  if( $test ) then
    echo "solvent scale converged"
    break
  endif

end

echo "capturing final solvent map in total_solvent${n}.map"
mapmask mapin total_solvent.map mapout total_solvent${n}.map << EOF >> $logfile
axis X Y Z
xyzlim asu
EOF
# this is the new master solvent
set mapfile = total_solvent${n}.map


set mtzfile = refmacout.mtz
set pdbfile = refmacout.pdb
converged:



# make sure map has right axes, etc.
mapmask mapin $mapfile mapout ${tempfile}mapfile.map << EOF >> $logfile
AXIS X Y Z
xyzlim asu
EOF





echo "collapsing $mtzfile supercell to ASU"
set reindex = `echo $super_mult | awk '{print "h/"$1",k/"$2",l/"$3}'`
reindex hklin $mtzfile hklout refmacout_small.mtz << EOF >! $logfile
reindex $reindex
symm $smallSGnum
EOF
echo "making $diffmapfile with extra fftB=$fftB"
fft hklin refmacout_small.mtz mapout ${tempfile}ffted.map <<EOF >> $logfile
scale F1 1 $fftB
labin F1=DELFWT PHI=PHDELWT
$GRID
EOF

set dens = `echo | mapdump mapin ${tempfile}ffted.map | awk '/density/{print $NF}'`
echo "$fftB $dens" | tee -a dens_vs_fftB.txt

# structure factor scales with number of unit cells
set cscale = `echo $ncells | awk '{print 1/$1}'`
mapmask mapin ${tempfile}ffted.map mapout $diffmapfile << EOF >> $logfile
scale factor $cscale
AXIS X Y Z
xyzlim asu
EOF



# ' "probability that something is there" '
# ' prob(rho) = sign(rho)*pow(erf(abs(rho/sigma(rho))/sqrt(2)),V/2/d^3) '
# ' add a scaled-down rho to break ties of zeroes? '

# number of shannon voxels expected
set exponent = `echo $smallCELLvolume $symops $solventB $fftB | awk '{V=$1;n=$2;B=$3+$4;pi=atan2(1,1)*4;print V/n/2/(sqrt((B+9.484)*log(2))/pi)**3}'`

echo "taking absolute value of Fo-Fc scaled to sigma=0.7071"
mapmask mapin $diffmapfile mapout fofc_sigma.map << EOF >> $logfile
scale sigma 0.70701 0
EOF
map_func.com -func abs  fofc_sigma.map -output abs.map >> $logfile
echo "erfpow"
map_func.com -func erfpow -param $exponent  abs.map -output significance.map >> $logfile

if( $maskB != 0 ) then
  echo "smoothing significance.map with scale=$mask_scale B=$maskB"
  map_scaleB_runme.com significance.map \
    scale=$mask_scale B=$maskB \
    reso=$mask_reso \
    output=${tempfile}offsetted.map >> $logfile

  echo "clip"
  map_func.com -func max -param 0 ${tempfile}offsetted.map -output ${tempfile}clip1.map >> $logfile
  map_func.com -func min -param 1 ${tempfile}clip1.map -output filtered_significance.map >> $logfile
  set significancemap = filtered_significance.map
else
  set significancemap = significance.map
endif


echo "making mask around protein in supercell"
make_supermask_runme.com $pdbfile $diffmapfile wetmask=0 maskB=50 outfile=protein_mask.map >> $logfile

echo "combining protein mask with significance mask into weight.map"
mapmask mapin1 $significancemap mapin2 protein_mask.map mapout weight.map << EOF >> $logfile
maps mult
EOF
echo "multiplying $diffmapfile by weight.map"
mapmask mapin1 weight.map mapin2 $diffmapfile mapout addme.map << EOF >> $logfile
maps mult
EOF

echo "adding weighted differences back to $mapfile "
mapmask mapin1 ${tempfile}mapfile.map mapin2 addme.map mapout posneg.map << EOF >> $logfile
maps add
EOF

echo "clipping off negative solvent density"
map_func.com -func max -param 0 posneg.map -output new_solvent.map >> $logfile

if(! $keepmaps) then
    rm -f total_solvent${n}.map start${n}_small.mtz >& /dev/null
endif
@ n = ( $n + 1 )
cp new_solvent.map total_solvent${n}.map
echo "solvent map is now total_solvent${n}.map"
set mapfile = total_solvent${n}.map

# measure how much map changed
map_func.com -func subtract $mapfile ${tempfile}mapfile.map \
   -output ${tempfile}diff_solvent.map >> $logfile
set stats = `echo | mapdump mapin ${tempfile}diff_solvent.map | awk '/dens/{print $NF}'`
echo "min max mean rms change in density: $stats"
set test = `echo $stats $min_delta_rho | awk '{d=sqrt($1*$1+$2*$2);print ( d < $NF )}'`
if( $test && $Bramp == 0 ) then
   echo "exiting because delta-rho < $min_delta_rho"
   break
endif

if( $skiprefine ) goto exit

set scaleoverall = 1
set Boverall = 0

end



exit:

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $?DEBUG && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif

if($?BAD) then
   echo "ERROR: $BAD"
   exit 9
endif
exit

#############################
## notes
#############################



Bramp:

set logfile = tempfile

cp scaled_solvent1.map 
set mapfile = total_solvent.map

cp scaled_solvent1.map $mapfile

foreach B ( `awk 'BEGIN{for(B=10000;B>2;B/=1.5)print sprintf("%.2g\n",B)+0}'` )

echo "B=$B"
touchup_solvent_runme.com refmacout.??? $mapfile skiprefine fftB=$B maskB=$B | tee touchup_B${B}.log

set mapfile = `awk '/^solvent mask is now/{print $NF}' touchup_B${B}.log`
solmap2refme.com $mapfile 

@ n = ( $n + 1 )
cp refmacout.pdb start${n}.pdb
echo "converge $n"
converge_refmac.com append nosalvage refme.mtz start${n}.pdb trials=5 >&! converge${n}.log

set refscale = `tail -n 1 refmac_scales.log | awk '{print $4+0}'`
set refB = `tail -n 1 refmac_scales.log | awk '{print $5+0}'`
map_scaleB_runme.com $mapfile scale=$refscale B=$refB output=total_solvent.map
set mapfile = total_solvent.map

end

exit

cd $pwd
foreach ramp ( `seq -f%02.0f 1 99` )
  cd $pwd
  mkdir -p ramp_${ramp}/touchup
  cd ramp_${ramp}/touchup
  cp ../noH/refmacout_minRfree.pdb .
  cp ../noH/scaled_solvent1.map .

  solmap2refme.com scaled_solvent1.map

  converge_refmac.com trials=1 refmacout_minRfree.pdb refme.mtz append nosalvage >&! converge1.log &

  set mapfile = `ls -1rt *solvent*.map | tail -n 1`
  echo $mapfile
  touchup_solvent_runme.com refmacout.??? $mapfile Bramp=1 >&! touchup_Bramp.log &

  sleep 1

  cd $pwd
end


foreach fftB ( `seq 0 1 10` `seq 15 5 50` `seq 60 10 100` `seq 200 100 1000` `seq 2000 1000 10000` )
fft hklin refmacout_small.mtz mapout ${tempfile}ffted.map <<EOF >! $logfile
scale F1 1 $fftB
labin F1=DELFWT PHI=PHDELWT
EOF
set dens = `awk '/density/{print $NF}' $logfile`
echo "$fftB $dens" | tee -a dens_vs_fftB.txt
end




if( 0 ) then
# rm  *
mapmask mapin ../rewater/total_solvent.map mapout total_solvent0.map << EOF > /dev/null
scale factor 0.020833333
axis X Y Z
xyzlim asu
EOF

cp ../rewater/salted_protein.pdb .
egrep -v "NH4|ACY|AMM|ACT" salted_protein.pdb >! protein.pdb

cp total_solvent0.map total_solvent.map
set scaleoverall = 0.85
set Boverall = 4

set s = 0
set n = 0
cp protein.pdb start0.pdb

touchup_solvent_runme.com start0.pdb total_solvent0.map \
  scale=0.85 B=4 iterate=100 mask_scale=1 maskB=5 >&! touchup_solvent_runme1.log &

endif

