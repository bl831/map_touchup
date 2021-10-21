#! /bin/tcsh -f
#
#   convert solvent map into partial structure for refmac
#   optionally allow subtraction of final, minimized solvent density
#       this is so that all coordinates can be refined
#
#
set pdbfile = ""
set mapfile = "solvent.map"
set outputmtz  = "refme.mtz"
set outputmap  = "Fpart.map"

set super_mult = ( 2 2 3 )
set reso = 0.95
set mapscale = 1
set proscale = 0
set negscale = 1

#set maskdown = "salt,water"
set subtract = solvent
set maskdown = none
set maskB = 5
set mask_reso = 1.3

set pdir = `dirname $0`
set path = ( $path $pdir )

set refmtz = ${pdir}/1aho.mtz

set tempfile = /dev/shm/${USER}/solmap_$$_
mkdir -p /dev/shm/${USER}

set logfile = solmap2refme_details.log

if(-e solmap_opts.txt) then
    echo "found solmap_opts.txt"
    cat solmap_opts.txt
    source solmap_opts.txt
endif


foreach arg ( $* )
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set Val = `echo $arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if("$Key" =~ *.pdb && "$val" == "") set pdbfile = $Key
    if("$key" == "pdbfile") set pdbfile = "$Val"
    if("$Key" =~ *.map && "$val" == "") set mapfile = $Key
    if("$key" == "mapfile") set mapfile = "$Val"
    if("$key" == "refmtz") set refmtz = "$Val"

    if("$key" == "outputmtz") set outputmtz = "$Val"
    if("$key" == "outputmap") set outputmap = "$Val"
    if("$key" == "reso") set reso = "$num"
    if("$key" == "mapscale") set mapscale = "$num"
    if("$key" == "proscale") set proscale = "$num"
    if("$key" == "negscale") set negscale = "$num"

    if("$key" == "subtract") set subtract = "$val"
    if("$key" == "maskdown") set maskdown = "$val"
    if("$Key" == "maskB") set maskB = "$num"
    if("$key" == "mask_reso") set mask_reso = "$num"

    if("$key" == "tempfile") then
        set tempfile = "$Val"
        setenv DEBUG
    endif
    if("$key" == "debug") then
        setenv DEBUG
    endif
end

set subtract = `echo $subtract | awk '{gsub("solvent","salt,water");print}'`
set maskdown = `echo $maskdown | awk '{gsub("solvent","salt,water");print}'`

set smallmtz = `echo $outputmtz | awk '{gsub(".mtz$","_small.mtz");print}'`

if(! -e "$mapfile") then
    set BAD = "cannot read $mapfile"
    goto exit
endif
if(! -e "$refmtz") then
    set BAD = "cannot read $refmtz"
    goto exit
endif
if(! -e "$pdbfile" && "$pdbfile" != "") then
    set BAD = "cannot read $pdbfile"
    goto exit
endif

cat << EOF
pdbfile= $pdbfile
mapfile= $mapfile

mapscale = $mapscale
proscale = $proscale
negscale = $negscale

subtract = $subtract
maskdown = $maskdown
maskB = $maskB
mask_reso = $mask_reso

EOF

# apply any scale factor and put in standard ASU
mapmask mapin $mapfile mapout ${tempfile}total_solvent0.map << EOF >> $logfile
scale factor $mapscale
xyzlim asu
axis X Y Z
EOF
set thismap = ${tempfile}total_solvent0.map

if(-e "$pdbfile") then
    set temp = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
    if($#temp == 6) set CELL = ( $temp )
endif

echo | mapdump mapin $thismap >! ${tempfile}mapdump.txt
set GRID = `awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}' ${tempfile}mapdump.txt`
set smallCELL = `awk '/Cell dimensions /{print $4,$5,$6,$7,$8,$9;exit}' ${tempfile}mapdump.txt`
set smallSGnum  = `awk '/ Space-group /{print $NF;exit}' ${tempfile}mapdump.txt`
if(! $?CELL) then
    set CELL = `echo $smallCELL $super_mult | awk '{print $1*$7,$2*$8,$3*$9,$4,$5,$6}'`
else
    set super_mult = `echo $CELL $smallCELL | awk '{printf("%.0f %.0f %.0f", $1/$7,$2/$8,$3/$9)}'`
endif
set ncells = `echo $super_mult | awk '{printf("%.0f\n",$1*$2*$3)}'`
set symops = `awk -v n=$smallSGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`

echo $smallCELL |\
awk 'NF==6{DTR=atan2(1,1)/45; A=cos(DTR*$4); B=cos(DTR*$5); G=cos(DTR*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf "%.3f\n", $1*$2*$3*sqrt(skew)}' |\
cat >! ${tempfile}volume
set smallCELLvolume = `cat ${tempfile}volume`
rm -f ${tempfile}volume
set CELLvolume = `echo $smallCELLvolume $ncells | awk '{print $1*$2}'`


set avgrho = `awk '/Mean density/{print $NF}' ${tempfile}mapdump.txt`
set bigwaters = `echo $avgrho 1 $CELL | awk '{print $1/10*$2*$3*$4*$5}'`
echo "equivalent of $bigwaters 10-electron supercell waters in ${mapscale}* $mapfile"

# split $pdbfile into protein and solvent
set negatoms = 0
set maskatoms = 0
set proatoms = 0
if(-e "$pdbfile") then
    # make PDBs of solvent only in big and small cells
    set pdbSG = `awk -v SGnum="$smallSGnum" -F "[\047]" '$1+0==SGnum{print $2;exit}' ${CLIBD}/symop.lib`
    echo $smallCELL $pdbSG |\
    awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' |\
    cat >! ${tempfile}smallcell.pdb
    echo $CELL P 1 |\
    awk '{printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s %s %s %s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)}' |\
    cat >! ${tempfile}bigcell.pdb

    # extract protein
    cp ${tempfile}smallcell.pdb ${tempfile}protein.pdb
    egrep "^ATOM|^HETAT" $pdbfile |\
    awk '{typ=substr($0,18,3)}\
      typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL/ || \
      typ~/HID|HIE|HIS|ILE|LEU|LYS|PHE|PRO|SER|THR|TRP|TYR/{\
       print}' >> ${tempfile}protein.pdb
    set proatoms = `egrep "^ATOM|^HETAT" ${tempfile}protein.pdb | wc -l`

    # extract salt
    cp ${tempfile}smallcell.pdb ${tempfile}salt.pdb
     egrep "^ATOM|^HETAT" $pdbfile |\
    awk '{typ=substr($0,18,3)}\
      typ~/ALA|ARG|ASN|ASP|ASH|CYS|CYX|GLN|GLU|GLH|GLY|VAL/{next}\
      typ~/HID|HIE|HIS|ILE|LEU|LYS|PHE|PRO|SER|THR|TRP|TYR|HOH/{next}\
       {print}' >> ${tempfile}salt.pdb
    set saltatoms = `egrep "^ATOM|^HETAT" ${tempfile}salt.pdb | wc -l`

     # extract waters
    cp ${tempfile}smallcell.pdb ${tempfile}water.pdb
    egrep "^ATOM|^HETAT" $pdbfile |\
    awk '{typ=substr($0,18,3)}\
      typ~/HOH|WAT/{\
       print}' >> ${tempfile}water.pdb
    set waters = `egrep "^ATOM|^HETAT" ${tempfile}water.pdb | wc -l`

    # gather atoms to subtract from the map
    if("$subtract" != "none" && "$negscale" != "0") then
        cp ${tempfile}smallcell.pdb ${tempfile}subtractme.pdb
        if( "$subtract"  =~ "*protein*") then
            echo "subtract: protein"
            egrep "^ATOM|^HETAT" ${tempfile}protein.pdb >> ${tempfile}subtractme.pdb
        endif
        if( "$subtract"  =~ "*salt*") then
            echo "subtract: salt"
            egrep "^ATOM|^HETAT" ${tempfile}salt.pdb >> ${tempfile}subtractme.pdb
        endif
        if( "$subtract"  =~ "*water*") then
            echo "subtract: water"
            egrep "^ATOM|^HETAT" ${tempfile}water.pdb >> ${tempfile}subtractme.pdb
        endif
        set negatoms = `egrep "^ATOM|^HETAT" ${tempfile}subtractme.pdb | wc -l`
   endif
 
    # gather atoms to squash in the map
    if("$maskdown" != "none") then
        cp ${tempfile}bigcell.pdb ${tempfile}maskme.pdb
        if( "$maskdown"  =~ "*protein*") then
            echo "maskdown: protein"
            egrep "^ATOM|^HETAT" ${tempfile}protein.pdb >> ${tempfile}maskme.pdb
        endif
        if( "$maskdown"  =~ "*salt*") then
            echo "maskdown: salt"
            egrep "^ATOM|^HETAT" ${tempfile}salt.pdb >> ${tempfile}maskme.pdb
        endif
        if( "$maskdown"  =~ "*water*") then
            echo "maskdown: water"
            egrep "^ATOM|^HETAT" ${tempfile}water.pdb >> ${tempfile}maskme.pdb
        endif
        set maskatoms = `egrep "^ATOM|^HETAT" ${tempfile}maskme.pdb | wc -l`
    endif
 endif


if( $negatoms && "$negscale" != "0") then
    echo "subtracting ${negscale}* $subtract in $pdbfile from ${mapscale}* $mapfile"

    sfall xyzin ${tempfile}subtractme.pdb mapout ${tempfile}sfalled.map << EOF >> $logfile
mode atmmap
CELL $smallCELL
SYMM $smallSGnum
SFSG 1
#BRESET 5
GRID $GRID
EOF
    set Fcscale = `echo $negscale $ncells $symops | awk '{print -$1/$2/$3}'`
    mapmask mapin ${tempfile}sfalled.map mapout ${tempfile}neg.map << EOF >> $logfile
scale factor $Fcscale
xyzlim asu
axis X Y Z
EOF
    set avgrho = `echo | mapdump mapin ${tempfile}neg.map | awk '/Mean density/{print -$NF}'`
    set bigwaters = `echo $avgrho 1 $CELL | awk '{print $1/8/$2*$3*$4*$5}'`
    echo "subtracting equivalent of $bigwaters waters from supercell"
    mapmask mapin1 $thismap mapin2 ${tempfile}neg.map \
       mapout ${tempfile}total_solvent1.map << EOF >> $logfile
maps add
EOF
    set thismap = ${tempfile}total_solvent1.map
endif


if( $proatoms && "$proscale" != "0") then
    echo "adding ${proscale}* protein in $pdbfile to the map"

    sfall xyzin ${tempfile}protein.pdb mapout ${tempfile}sfalled.map << EOF >> $logfile
mode atmmap
CELL $smallCELL
SYMM $smallSGnum
SFSG 1
#BRESET 5
GRID $GRID
EOF
    set Fcscale = `echo $proscale $ncells $symops | awk '{print $1/$2/$3}'`
    mapmask mapin ${tempfile}sfalled.map mapout ${tempfile}pro.map << EOF >> $logfile
scale factor $Fcscale
xyzlim asu
axis X Y Z
EOF
    mapmask mapin1 $thismap mapin2 ${tempfile}pro.map \
       mapout ${tempfile}total_density.map << EOF >> $logfile
maps add
EOF
    set thismap = ${tempfile}total_density.map
endif



if("$maskdown" != "none" && -e ${tempfile}maskme.pdb ) then
    echo "squashing density around $maskdown from $pdbfile with mask"
    make_supermask_runme.com ${tempfile}maskme.pdb $thismap maskB=$maskB \
     wetmask=1 outfile=${tempfile}smallmask.map >> $logfile
    mapmask mapin1 $thismap \
      mapin2 ${tempfile}smallmask.map \
      mapout ${tempfile}maskedsol.map << EOF >> $logfile
maps mult
EOF
    set thismap = ${tempfile}maskedsol.map
endif


if(-e ${tempfile}maskedsol.map) cp ${tempfile}maskedsol.map maskedsol.map
if(-e ${tempfile}smallmask.map) cp ${tempfile}smallmask.map mask.map
if(-e ${tempfile}neg.map) then
    echo scale factor -1 |\
      mapmask mapin ${tempfile}neg.map mapout pdbsolvent.map >> /dev/null
endif
cp $thismap $outputmap

set avgrho = `echo | mapdump mapin $outputmap | awk '/Mean density/{print $NF}'`
set bigwaters = `echo $avgrho 1 $CELL | awk '{print $1/10/$2*$3*$4*$5}'`
echo "equivalent of $bigwaters 10-electron supercell waters in  $outputmap"

echo "turning $outputmap into Fpart"
mapmask mapin $outputmap mapout ${tempfile}sfallme.map << EOF >> $logfile
xyzlim cell
axis Z X Y
EOF
sfall mapin ${tempfile}sfallme.map hklout ${tempfile}sfalled.mtz << EOF >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF

echo "adding Fpart to $refmtz"
cad hklin1 $refmtz hklin2 ${tempfile}sfalled.mtz hklout $smallmtz << EOF >> $logfile
labin file 1 all
labin file 2 E1=FC E2=PHIC
labou file 2 E1=Fpart E2=PHIpart
EOF

echo "expanding to supercell in $outputmtz"
cad hklin1 $smallmtz hklout ${tempfile}.mtz << EOF >> $logfile
labin file 1 all
outlim space 1
EOF
cad hklin1 ${tempfile}.mtz hklout ${tempfile}P1.mtz << EOF  >> $logfile
labin file 1 all
scale file 1 $ncells 0
symm 1
EOF
echo $super_mult |\
awk '{print "reindex "$1"h,"$2"k,"$3"l"}' |\
 reindex hklin ${tempfile}P1.mtz hklout $outputmtz >> $logfile


#ls -l $outputmap
#ls -l $smallmtz
#ls -l $outputmtz


exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if(! $?DEBUG && "$tempfile" != "") then
    echo "cleaning up..."
    rm -rf ${tempfile}* >& /dev/null
endif

exit

###############################################################
###############################################################
###############################################################


# check that phases really do work after re-indexing


mapmask mapin solvent.map mapout sfallme.map << EOF
xyzlim cell
axis Z X Y
EOF

sfall mapin sfallme.map hklout sfalled.mtz << EOF
mode sfcalc mapin
resolution $reso
sfsg 1
EOF


echo reindex h/2,k/2,l/3 |\
 reindex hklin erefmac_Prod/refme.mtz hklout temp.mtz
cad hklin1 temp.mtz hklout this.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthis E2=PHIthis
symm P212121
EOF
cad hklin1 refmac_Prod/refme.mtz hklout that.mtz << EOF
labin file 1 E1=Fpart E2=PHIpart
labou file 1 E1=Fthat E2=PHIthat
symm P212121
EOF
cad hklin1 this.mtz hklin2 sfalled.mtz hklout phases.mtz << EOF
labin file 1 E1=PHIthis
labin file 2 E1=PHIC
EOF
mtz2txt phases.mtz 

awk '$4~/[0-9]/ && $5~/[0-9]/{print sqrt(($4-$5)^2),$0}' phases.csh  | sort -g




