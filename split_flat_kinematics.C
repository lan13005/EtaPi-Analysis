bool filterOmega(float omega, float Mpi0eta){
    // omega should be in degrees and mass in GeV
    return -29.0*atan(-1.05*Mpi0eta+2.78)+328 > omega;
}

void split_flat_kinematics(){
    // ********************************************
    // THIS PROGRAM WILL SPLIT ROOT FILES INTO DIFFERENT BEAM ENERGY + T BINNINGS
    //    BOTH THE THROWN AND RECON VALUES WILL BE USED IN DETERMINING BINNING 
    //    IF THEY EXISTS. 
    //    I.E. RECONSTRUCTED MC WOULD HAVE BOTH.
    //         DATA WOULD ONLY HAVE RECON
    //         THROWN MC WOULD HAVE ONLY THROWN VALUES
    //    BY DEFAULT THERE IS ALWAYS A BINNING INTO SEPARATE POLARIZATIONS WHICH
    //         USED FOR AMPLITUDE FITS WITH POLARIZATION
    // ********************************************

    bool sumRuns=true;
    bool forceSplitting=true; // Should we run the splitting again? Or should we just sum runs if sumRuns=true
    bool remergePols=true; // should we remerge polarizations after splitting? 
    
    //string extraTag="_mpipGT20_selectGenT";
    string extraTag="";

    string folder="phase1_selected_v2/";
    bool ignorePolarization=false;
    vector<string> runs={"2017_1","2018_1","2018_8"};
    vector<string> files;
    for (auto run: runs){
        //files.push_back("D"+run+"_selected_data_flat.root");
        //files.push_back("D"+run+"_selected_bkgnd_flat.root");
        //files.push_back("D"+run+"_selected_acc_flat.root");
        //files.push_back("F"+run+"_selected_acc_flat.root");
        files.push_back("F"+run+"_gen_data_flat.root");
    }

//    string folder="kmatrix_selected_v1/";
//    bool ignorePolarization=true; // if true then beamAngle will be set to 0
//    vector<string> runs={""};
//    vector<string> files;
//    for (auto run: runs){
//        files.push_back("kmatrix_selected_data_flat.root");
//        files.push_back("kmatrix_selected_bkgnd_flat.root");
//        files.push_back("kmatrix_gen_data_flat.root");
//        files.push_back("F2018_8_selected_acc_flat.root");
//        files.push_back("F2018_8_gen_data_flat.root");
//    }

    // ********************************************
    // HERE IS ALL THE BINNINGS WE WILL LOOP OVER
    // ********************************************
    int nFileTypes=((int)files.size())/((int)runs.size());

    map<int,int> pols={{0,0},{45,1},{90,2},{135,3},{-1,4}};
    vector<string> polstrings={"000","045","090","135","AMO"};

    map<string,int> ts={{"010020",0},{"0200325",1},{"0325050",2},{"050075",3},{"075100",4}}; // t
    vector<float> mint={0.1,0.2,0.325,0.5,0.75};
    vector<float> maxt={0.2,0.325,0.5,0.75,1.0};
    map<string,int> mpi0etas={{"104180",0}}; // m 
    vector<float> minmpi0eta={1.04};
    vector<float> maxmpi0eta={1.80};
    //map<string,int> ts={{"all",0}}; // t
    //vector<float> mint={0};
    //vector<float> maxt={100};
    //map<string,int> mpi0etas={{"104156",0}}; // m 
    //vector<float> minmpi0eta={1.04};
    //vector<float> maxmpi0eta={1.56};

    const int nts=(const int)mint.size();
    const int nms=(const int)minmpi0eta.size();

    // ********************************************
    // CREATE SOME EXTRA FOLDERS FOR ORGANIZATION
    // ********************************************
    
    if (forceSplitting){
        for (auto const& t: ts){
            for (auto const& m: mpi0etas){
                string floc=folder+"t"+t.first+"_m"+m.first+extraTag+"/";
                gSystem->Exec(("mkdir -p "+floc).c_str()); } }

        for (auto file: files){
            cout << "================================================" << endl;
            cout << "FILE: " << file << endl;
            cout << "================================================" << endl;

            TFile *oldfile = new TFile((folder+file).c_str());
            TTree *oldtree = (TTree*)oldfile->Get("kin");
            Long64_t nentries = oldtree->GetEntries();

            // ********************************************
            /////// CREATE A BUNCH OF ROOT FILES AND TREES
            // ********************************************
            // The order of the array dimensions represents: #pols, #t, #e, #m
            TFile *newfile[5][nts][nms];
            TTree *newtree[5][nts][nms]; 
            int ip, it, im;
            ip=0;
            for (auto const& pol: pols){ it=0;
                for (auto const& t: ts){ im=0;
                    for (auto const& m: mpi0etas){
                        string floc=folder+"t"+t.first+"_m"+m.first+extraTag+"/";
                        newfile[ip][it][im] = new TFile((floc+"pol"+polstrings[ip]+"_t"+t.first+"_m"+m.first+extraTag+"_"+file).c_str(),"recreate");
                        newtree[ip][it][im] = oldtree->CloneTree(0);
                        ++im;
                    } ++it;
                } ++ip;
            }

            // ********************************************
            /////// CHECKING TO WHICH BRANCHES EXIST - WHETHER DATA OR SIMULATIONS
            // ********************************************
            int BeamAngle;
            int beamAngle=0; // this is the angle that is actually used but will depend on ignorePolarization variable
            float mandelstam_t;
            float Ebeam;
            float mpi0eta;
            float mpi0p;
            float mandelstam_t_thrown;
            float Ebeam_thrown;
            float mpi0eta_thrown;
            float vanHove_omega;
            bool pVH_pi0p;
            bool pVH_pi0p2;
            // Check to see if we have the right branches for the different data,bkgnd,acc,gen trees
            //    Small caveat - thrown branches ALWAYS exist even for data,bkgnd trees but the branch will only contain 0s
            bool has_recon_branches = (bool)oldtree->GetListOfBranches()->FindObject("Ebeam"); // returned object decays to a boolean
            cout << "has proper recon branches: " << has_recon_branches << endl;
            
            oldtree->SetBranchAddress("BeamAngle",&BeamAngle);
            if (has_recon_branches){
                oldtree->SetBranchAddress("mandelstam_t",&mandelstam_t);
                oldtree->SetBranchAddress("Ebeam",&Ebeam);
                oldtree->SetBranchAddress("Mpi0eta",&mpi0eta);
                oldtree->SetBranchAddress("Mpi0p",&mpi0p);
                oldtree->SetBranchAddress("pVH_pi0p",&pVH_pi0p);
                oldtree->SetBranchAddress("vanHove_omega",&vanHove_omega);
            }
            oldtree->SetBranchAddress("mandelstam_t_thrown",&mandelstam_t_thrown);
            oldtree->SetBranchAddress("Ebeam_thrown",&Ebeam_thrown);
            oldtree->SetBranchAddress("Mpi0eta_thrown",&mpi0eta_thrown);
            oldtree->GetEntry(0);
            bool has_thrown_branches = (Ebeam_thrown!=0);
            cout << "has proper thrown branches: " << has_thrown_branches << endl;

            bool isData=!has_thrown_branches*has_recon_branches; // does not have thrown but has recon branches
            bool isAcc=has_recon_branches*has_thrown_branches; // has both recon and thrown branches
            bool isGen=!has_recon_branches*has_thrown_branches; // does not have recon branches but has thrown branches

            // ********************************************
            /////// FILL TREES IN THE CORRECT BINNING
            // ********************************************
            for (Long64_t i=0;i<nentries; i++) {
                 oldtree->GetEntry(i);
                 it=0;
                 if (!ignorePolarization)
                     beamAngle=BeamAngle;
                 for (auto const& t: ts){ im=0;
                     for (auto const& m: mpi0etas){
                         //pVH_pi0p2=filterOmega(vanHove_omega,mpi0eta);
                         //if (has_recon_branches*!pVH_pi0p2) continue;
                         //if (has_recon_branches*!(mpi0p>2.0)) continue;
                         if (has_recon_branches*!((mandelstam_t>mint[it])*(mandelstam_t<maxt[it]))) continue;
                         if (has_recon_branches*!((mpi0eta>minmpi0eta[im])*(mpi0eta<maxmpi0eta[im]))) continue;
                         //if (has_thrown_branches*!((mandelstam_t_thrown>mint[it])*(mandelstam_t_thrown<maxt[it]))) continue;
                         //if (has_thrown_branches*!((mpi0eta_thrown>minmpi0eta[im])*(mpi0eta_thrown<maxmpi0eta[im]))) continue;
                         //if (isGen*!((mandelstam_t_thrown>mint[it])*(mandelstam_t_thrown<maxt[it]))) continue;
                          newtree[pols[beamAngle]][it][im]->Fill();
                          ++im;
                     } ++it;
                 }
            }


            // ********************************************
            /////// WRITE TREES TO THEIR ASSOCIATED FILES
            // ********************************************
            int post_nentries=0;
            int tmp_nentries;
            for (int ip=0; ip<(int)pols.size(); ++ip){
                for (int it=0; it<(int)ts.size(); ++it){
                    for (int im=0; im<(int)mpi0etas.size(); ++im){
                        newfile[ip][it][im]->cd();
                        newtree[ip][it][im]->Write();
                        tmp_nentries=newtree[ip][it][im]->GetEntries();
                        post_nentries+=tmp_nentries;
                        cout << "Entries in [pol,tbin,mbin=" << polstrings[ip] << "," << it << "," << im << "] = " << tmp_nentries << endl;
                        newfile[ip][it][im]->Close();
            }}}
        } // close the forceSplitting condition
    } // closes files loop


    // ********************************************
    /////// COMPLICATED FUNCTION TO HADD FILES FROM DIFFERENT RUNS
    // We will loop through all runs and create new files with TOT replacing the run tag
    // We do this by creating a giant hadd string and executing the system command
    // ********************************************
    int ip, it, im;
    if (sumRuns){
        it=0;
        for (auto const& t: ts){ im=0;
            for (auto const& m: mpi0etas){ 
                for (int j=0; j<nFileTypes; ++j){ ip=0;
                    string remergePolCmd="hadd -f ";
                    remergePolCmd+=folder+"t"+t.first+"_m"+m.first+extraTag+"/"+"polALL_t"+t.first+"_m"+m.first+extraTag+"_";
                    remergePolCmd+=files[j][0]+(string)"TOT"+files[j].substr(runs[0].size()+1,files[j].size());
                    for (auto const& pol: pols){
                        string cmd;
                        string target;
                        for (int i=0; (int)i<runs.size(); ++i){
                            if (i==0){
                                target=folder+"t"+t.first+"_m"+m.first+extraTag+"/"+"pol"+polstrings[ip]+"_t"+t.first+"_m"+m.first+extraTag+"_";
                                target+=files[j][0]+(string)"TOT"+files[j].substr(runs[0].size()+1,files[j].size());
                                cmd="hadd -f "+target;
                            }
                            cmd+=" "+folder+"t"+t.first+"_m"+m.first+extraTag+"/"+"pol"+polstrings[ip]+"_t"+t.first+"_m"+m.first+extraTag+"_"+files[nFileTypes*i+j];
                        }
                        cout << endl << cmd << endl;
                        gSystem->Exec(cmd.c_str()); 
                        remergePolCmd+=" "+target;
                        ++ip;
                    }
                    if (remergePols){
                        cout << "\n---- REMERGING POLS ----\n" << remergePolCmd << endl;
                        gSystem->Exec(remergePolCmd.c_str()); 
                        cout << endl;
                    }
                } 
                ++im;
            } ++it;
        }
    }
}








