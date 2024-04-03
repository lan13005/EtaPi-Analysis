#!/usr/bin/perl

use Cwd;

$lowMass = 1.04;
$highMass = 1.56;
$nBins =  13;
$fitName = "EtaPi_fit";

$maxEvts = 1E9;

$workingDir=getcwd();
print "\n\ncurrent working dir: $workingDir";
print "\n===================================\n";

$t="010020";
$baseGenDir="/w/halld-scshelf2101/lng/WORK/EtaPi-Analysis/study_pwa/mass_independent_fits/rootFiles/t$t\_m104180_selectGenTandM/";
$baseAccDir="/w/halld-scshelf2101/lng/WORK/EtaPi-Analysis/study_pwa/mass_independent_fits/rootFiles/t$t\_m104172_selectGenTandM_nominal/";
$baseBkgDir="/w/halld-scshelf2101/lng/WORK/EtaPi-Analysis/study_pwa/mass_independent_fits/rootFiles/t$t\_m104172_selectGenTandM_nominal/";
$baseDatDir="/w/halld-scshelf2101/lng/WORK/EtaPi-Analysis/study_pwa/mass_independent_fits/rootFiles/t$t\_m104172_selectGenTandM_nominal/";
$baseGenFileName="_t$t\_m104180_selectGenTandM_FTOT_gen_data_flat"; # cannot contain the file extension .root
$baseAccFileName="_t$t\_m104172_selectGenTandM_FTOT_selected_nominal_acc_flat";
$baseBkgFileName="_t$t\_m104172_selectGenTandM_DTOT_selected_nominal_bkgnd_flat";
$baseDatFileName="_t$t\_m104172_selectGenTandM_DTOT_selected_nominal_data_flat";

@polTags_DB=qw(000 045 090 135); # D/B = data, bkgnd
@polTags_AG=qw(000 045 090 135); # A/G = accmc, genmc. i.e. we can choose qw(ALL)

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
mkdir $fitDir unless -d $fitDir;
#### CHOICE 2 : Create a ram disk if you have access to lots of ram
#system( "rm -f $fitName" );
#system( "rm -rf /dev/shm/$fitName" );
#system( "mkdir /dev/shm/$fitName" );
#system( "ln -s /dev/shm/$fitName $fitName" );

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

