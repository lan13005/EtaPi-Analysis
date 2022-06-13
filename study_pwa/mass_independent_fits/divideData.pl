#!/usr/bin/perl

use Cwd;

$lowMass = 0.80;#1.04;
$highMass = 1.80;
$nBins =  25;#19;
$fitName = "EtaPi_fit";

# put a limit on the number of data events to process
# gen MC and acc MC smaples are not limited
$maxEvts = 1E9;

$workingDir=getcwd();
print "\n\ncurrent working dir: $workingDir";
print "\n===================================\n";

$t="0103";
$m="104180";
#$baseGenDir="/d/grid17/ln16/dselector_v2/test/phase1_selected/t$t\_m$m/";
#$baseAccDir="/d/grid17/ln16/dselector_v2/test/phase1_selected/t$t\_m$m/";
#$baseBkgDir="/d/grid17/ln16/dselector_v2/test/phase1_selected/t$t\_m$m/";
#$baseDatDir="/d/grid17/ln16/dselector_v2/test/phase1_selected/t$t\_m$m/";
#$baseGenFileName="_t0103_m$m\_FTOT_gen_data_flat"; # cannot contain the file extension .root
#$baseAccFileName="_t0103_m$m\_FTOT_selected_acc_flat";
#$baseBkgFileName="_t0103_m$m\_DTOT_selected_bkgnd_flat";
#$baseDatFileName="_t0103_m$m\_DTOT_selected_data_flat";
$baseAccDir="/d/grid17/ln16/dselector_v2/test/kmatrix_selected/tall_m080180/";
$baseGenDir="/d/grid17/ln16/dselector_v2/test/kmatrix_selected/tall_m080180/";
$baseDatDir="/d/grid17/ln16/dselector_v2/test/kmatrix_selected/";
$baseBkgDir="/d/grid17/ln16/dselector_v2/test/kmatrix_selected/";
$baseAccFileName="_tall_m080180_F2018_8_selected_acc_flat";
$baseGenFileName="_tall_m080180_F2018_8_gen_data_flat";
$baseDatFileName="_tall_m080180_kmatrix_selected_halved_data_flat";
$baseBkgFileName="_tall_m080180_kmatrix_selected_halved_bkgnd_flat";

@polTags=qw(000);# 045 090 135);
print "DATAFILES:\n";
foreach $polTag (@polTags){
    print "$baseDatDir\pol$polTag$baseDatFileName\.root\n";
}
print "------------------\n";

print "BKGNDFILES:\n";
foreach $polTag (@polTags){
    print "$baseBkgDir\pol$polTag$baseBkgFileName\.root\n";
}
print "------------------\n";

print "ACCFILES:\n";
foreach $polTag (@polTags){
    print "$baseAccDir\pol$polTag$baseAccFileName\.root\n";
}
print "------------------\n";

print "GENFILES:\n";
foreach $polTag (@polTags){
    print "$baseGenDir\pol$polTag$baseGenFileName\.root\n";
}
print "------------------\n";


$cfgTempl = "$workingDir/config_files/zlm_etapi_bothReflect_bothM.cfg";

### things below here probably don't need to be modified

# this is where the goodies for the fit will end up
$fitDir = "$workingDir/$fitName/";
print "Output fitDir: $fitDir";
print "\n";

#### CHOICE 1 : Create a regular folder 
#mkdir $fitDir unless -d $fitDir;
#### CHOICE 2 : Create a ram disk if you have access to lots of ram
system( "rm -f $fitName" );
system( "rm -rf /dev/shm/$fitName" );
system( "mkdir /dev/shm/$fitName" );
system( "ln -s /dev/shm/$fitName $fitName" );

chdir $fitDir;

print "Changing into $fitDir\n";

# use the split_mass command line tool to divide up the
foreach $polTag (@polTags){
    $fileTag="pol$polTag$baseDatFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseDatDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseDatDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );

    $fileTag="pol$polTag$baseBkgFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseBkgDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseBkgDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );

    $fileTag="pol$polTag$baseAccFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseAccDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseAccDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );

    $fileTag="pol$polTag$baseGenFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseGenDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseGenDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );
}

# make directories to perform the fits in
for( $i = 0; $i < $nBins; ++$i ){

  mkdir "bin_$i" unless -d "bin_$i";
  
  system( "mv *\_$i.root bin_$i" );

  chdir "bin_$i";

#we are essentially copying fit_etapi_moments.cfg and substituting some variables. CFGOUT is going to be a config file in all of our bins. CFGIN is fit_etapi_moments.cfg. Note how fit_etapi_moments.cfg has these place holders defined (DATAFILE,ACCMCFILE,GENMCFILE ... ). They will get replaced here to fit the bin directory. 
  open( CFGOUT, ">bin_$i-full.cfg" );
  open( CFGIN, $cfgTempl ); 

  while( <CFGIN> ){
    foreach $polTag (@polTags){
        s:DATAFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseDatFileName\_$i.root:;
        s:BKGNDFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseBkgFileName\_$i.root:;
        s:ACCMCFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseAccFileName\_$i.root:;
        s:GENMCFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseGenFileName\_$i.root:;
        s:NIFILE_$polTag:bin_$i\_$polTag.ni:;
    }

    s/FITNAME/bin_$i/;

    print CFGOUT $_;
  }

  close CFGOUT;
  close CFGIN;
  
  #system( "touch param_init.cfg" );

  chdir $fitDir;
}

