
# helper script to generate a small sample of DUNE 10kt detsim events for testing
N_EVENTS=10
lar -c prod_eminus_1-500MeV_dune10kt_1x2x6.fcl -o gen_eminus_1-500MeV_dune10kt_1x2x6.root -n ${N_EVENTS}
lar -c standard_g4_dune10kt_1x2x6.fcl -o g4_eminus_1-500MeV_dune10kt_1x2x6.root gen_eminus_1-500MeV_dune10kt_1x2x6.root -n -1
lar -c detsim_dune10kt_1x2x6_notpcsigproc.fcl -o detsim_eminus_1-500MeV_dune10kt_1x2x6.root g4_eminus_1-500MeV_dune10kt_1x2x6.root -n -1


# Run the trigger simulation and TTree writing
# lar -c triggersim_1x2x6_simpleThr_simpleWin_simpleWin.fcl -o trgsim_eminus_1-500MeV_dune10kt_1x2x6.root detsim_eminus_1-500MeV_dune10kt_1x2x6.root 
# lar -c triggerana_tree_1x2x6.fcl -T anatree_eminus_1-500MeV_dune10kt_1x2x6.root trgsim_eminus_1-500MeV_dune10kt_1x2x6.root