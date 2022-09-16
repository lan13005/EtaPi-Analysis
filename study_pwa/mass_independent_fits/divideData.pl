#!/usr/bin/perl

use Cwd;

#$genmc_filter=" mandelstam_t_thrown 0.1 1.0";
#$accmc_filter="";
#$bkgnd_filter="";
#$data_filter="";

$lowMass = 1.04;
$highMass = 1.80;
$nBins =  19;
$fitName = "EtaPi_fit";

# put a limit on the number of data events to process
# gen MC and acc MC smaples are not limited
$maxEvts = 1E9;

$workingDir=getcwd();
print "\n\ncurrent working dir: $workingDir";
print "\n===================================\n";

$t="010020";
$m="104180";
$extraTag="_vh_selectGenT";
$baseGenDir="/d/grid17/ln16/dselector_v3/phase1_selected_v2/t$t\_m$m$extraTag/";
$baseAccDir="/d/grid17/ln16/dselector_v3/phase1_selected_v2/t$t\_m$m$extraTag/";
$baseBkgDir="/d/grid17/ln16/dselector_v3/phase1_selected_v2/t$t\_m$m$extraTag/";
$baseDatDir="/d/grid17/ln16/dselector_v3/phase1_selected_v2/t$t\_m$m$extraTag/";
$baseGenFileName="_t$t\_m$m$extraTag\_FTOT_gen_data_flat"; # cannot contain the file extension .root
$baseAccFileName="_t$t\_m$m$extraTag\_FTOT_selected_acc_flat";
$baseBkgFileName="_t$t\_m$m$extraTag\_DTOT_selected_bkgnd_flat";
$baseDatFileName="_t$t\_m$m$extraTag\_DTOT_selected_data_flat";

#$baseAccDir="/d/grid17/ln16/dselector_v3/kmatrix_selected_v1/tall_m104156/";
#$baseGenDir="/d/grid17/ln16/dselector_v3/kmatrix_selected_v1/tall_m104156/";
#$baseDatDir="/d/grid17/ln16/dselector_v3/kmatrix_selected_v1/tall_m104156/";
#$baseBkgDir="/d/grid17/ln16/dselector_v3/kmatrix_selected_v1/tall_m104156/";
#$baseAccFileName="_tall_m104156_F2018_8_selected_acc_flat";
#$baseGenFileName="_tall_m104156_F2018_8_gen_data_flat";
#$baseDatFileName="_tall_m104156_kmatrix_selected_data_flat";
#$baseBkgFileName="_tall_m104156_kmatrix_selected_bkgnd_flat";

@polTags_DB=qw(000 045 090 135); # D/B = data, bkgnd
@polTags_AG=qw(000 045 090 135); # A/G = accmc, genmc. i.e. we can choose qw(ALL)
#@polTags_DB=qw(000); # D/B = data, bkgnd
#@polTags_AG=qw(000); # A/G = accmc, genmc. i.e. we can choose qw(ALL)

print "DATAFILES:\n";
foreach $polTag (@polTags_DB){
    print "$baseDatDir\pol$polTag$baseDatFileName\.root\n";
}
print "------------------\n";

print "BKGNDFILES:\n";
foreach $polTag (@polTags_DB){
    print "$baseBkgDir\pol$polTag$baseBkgFileName\.root\n";
}
print "------------------\n";

print "ACCFILES:\n";
foreach $polTag (@polTags_AG){
    print "$baseAccDir\pol$polTag$baseAccFileName\.root\n";
}
print "------------------\n";

print "GENFILES:\n";
foreach $polTag (@polTags_AG){
    print "$baseGenDir\pol$polTag$baseGenFileName\.root\n";
}
print "------------------\n";


$cfgTempl = "$workingDir/config_files/zlm_etapi_bothReflect_bothM_loop.cfg";

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
foreach $polTag (@polTags_DB){
    $fileTag="pol$polTag$baseDatFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseDatDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseDatDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );

    $fileTag="pol$polTag$baseBkgFileName";
    $dataFile="$fileTag\.root";
    print( "split_mass $baseBkgDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin\n" );
    system( "split_mass $baseBkgDir$dataFile $fileTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );
}

foreach $polTag (@polTags_AG){
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
    #s:ROOTDataReaderFilter\ LOOPDATA:ROOTDataReaderFilter\ LOOPDATA$data_filter:;
    #s:ROOTDataReaderFilter\ LOOPBKGND:ROOTDataReaderFilter\ LOOPBKGND$bkgnd_filter:;
    #s:ROOTDataReaderFilter\ LOOPGENMC:ROOTDataReaderFilter\ LOOPGENMC$genmc_filter:;
    #s:ROOTDataReaderFilter\ LOOPACCMC:ROOTDataReaderFilter\ LOOPACCMC$accmc_filter:;
    foreach $polTag (@polTags_DB){
        s:DATAFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseDatFileName\_$i.root:;
        s:BKGNDFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseBkgFileName\_$i.root:;
        if ($polTags_AG[0] eq "ALL"){
            s:ACCMCFILE_$polTag:${fitDir}bin_$i/polALL$baseAccFileName\_$i.root:;
            s:GENMCFILE_$polTag:${fitDir}bin_$i/polALL$baseGenFileName\_$i.root:;
        } else {
            s:ACCMCFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseAccFileName\_$i.root:;
            s:GENMCFILE_$polTag:${fitDir}bin_$i/pol$polTag$baseGenFileName\_$i.root:;
        }
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

