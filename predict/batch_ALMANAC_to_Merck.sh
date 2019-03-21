#!/bin/sh

# When running this script via a scheduler, you will need to specify the starting directory.
# To do this, you can add an explict "cd <starting directory>" here. This starting directory is
# assumed to be the directory that contains the "predict_synergy" program.

DRUG_FEATURES=../preprocess/ALMANAC_and_Merck.obable.FP2.csv
CELL_FEATURES=../preprocess/CCLE_DepMap_18q3_RNAseq_RPKM.LINCS1000.tsv

DATA_DIR=../preprocess/data_bliss

# Minimum synergy threshold
SYNERGY_THRESHOLD=-0.4

# The number of trees to use in the random forest algorithm
FOREST_SIZE=10000

# Use four ranks per node. Note that this assumes a Torque scheduler. If you are using a 
# different scheduler, you will need to modify how the number of nodes (i.e. the number
# of MPI ranks) is computed.
NUM_NODES=$(( $PBS_NUM_NODES * 4 ))

# The number of threads per rank to use.
export OMP_NUM_THREADS=12

mpirun -np $NUM_NODES \
	--oversubscribe \
	--bind-to socket \
	-x OMP_NUM_THREADS \
	./predict_synergy \
	--forest.size $FOREST_SIZE \
	--cell $CELL_FEATURES \
	--drug $DRUG_FEATURES \
	--synergy.test.random \
	--enable-pair-prediction \
	-o transfer.cell.drug.FP2.shared.$FOREST_SIZE.out \
	--synergy.threshold $SYNERGY_THRESHOLD \
	--synergy.train "CCLE.786O=$DATA_DIR/ALMANAC.CCLE_786O.csv" \
	--synergy.train "CCLE.A498=$DATA_DIR/ALMANAC.CCLE_A498.csv" \
	--synergy.train "CCLE.A549=$DATA_DIR/ALMANAC.CCLE_A549.csv" \
	--synergy.train "CCLE.ACHN=$DATA_DIR/ALMANAC.CCLE_ACHN.csv" \
	--synergy.train "CCLE.BT549=$DATA_DIR/ALMANAC.CCLE_BT549.csv" \
	--synergy.train "CCLE.CAKI1=$DATA_DIR/ALMANAC.CCLE_CAKI1.csv" \
	--synergy.train "CCLE.DU145=$DATA_DIR/ALMANAC.CCLE_DU145.csv" \
	--synergy.train "CCLE.EKVX=$DATA_DIR/ALMANAC.CCLE_EKVX.csv" \
	--synergy.train "CCLE.HCT116=$DATA_DIR/ALMANAC.CCLE_HCT116.csv" \
	--synergy.train "CCLE.HCT15=$DATA_DIR/ALMANAC.CCLE_HCT15.csv" \
	--synergy.train "CCLE.HL60=$DATA_DIR/ALMANAC.CCLE_HL60.csv" \
	--synergy.train "CCLE.HOP62=$DATA_DIR/ALMANAC.CCLE_HOP62.csv" \
	--synergy.train "CCLE.HOP92=$DATA_DIR/ALMANAC.CCLE_HOP92.csv" \
	--synergy.train "CCLE.HS578T=$DATA_DIR/ALMANAC.CCLE_HS578T.csv" \
	--synergy.train "CCLE.HT29=$DATA_DIR/ALMANAC.CCLE_HT29.csv" \
	--synergy.train "CCLE.IGROV1=$DATA_DIR/ALMANAC.CCLE_IGROV1.csv" \
	--synergy.train "CCLE.K562=$DATA_DIR/ALMANAC.CCLE_K562.csv" \
	--synergy.train "CCLE.KM12=$DATA_DIR/ALMANAC.CCLE_KM12.csv" \
	--synergy.train "CCLE.LOXIMVI=$DATA_DIR/ALMANAC.CCLE_LOXIMVI.csv" \
	--synergy.train "CCLE.MALME3M=$DATA_DIR/ALMANAC.CCLE_MALME3M.csv" \
	--synergy.train "CCLE.MCF7=$DATA_DIR/ALMANAC.CCLE_MCF7.csv" \
	--synergy.train "CCLE.MDAMB231=$DATA_DIR/ALMANAC.CCLE_MDAMB231.csv" \
	--synergy.train "CCLE.MDAMB435S=$DATA_DIR/ALMANAC.CCLE_MDAMB435S.csv" \
	--synergy.train "CCLE.MDAMB468=$DATA_DIR/ALMANAC.CCLE_MDAMB468.csv" \
	--synergy.train "CCLE.NCIH226=$DATA_DIR/ALMANAC.CCLE_NCIH226.csv" \
	--synergy.train "CCLE.NCIH23=$DATA_DIR/ALMANAC.CCLE_NCIH23.csv" \
	--synergy.train "CCLE.NCIH460=$DATA_DIR/ALMANAC.CCLE_NCIH460.csv" \
	--synergy.train "CCLE.NCIH522=$DATA_DIR/ALMANAC.CCLE_NCIH522.csv" \
	--synergy.train "CCLE.NIHOVCAR3=$DATA_DIR/ALMANAC.CCLE_NIHOVCAR3.csv" \
	--synergy.train "CCLE.OVCAR4=$DATA_DIR/ALMANAC.CCLE_OVCAR4.csv" \
	--synergy.train "CCLE.OVCAR8=$DATA_DIR/ALMANAC.CCLE_OVCAR8.csv" \
	--synergy.train "CCLE.PC3=$DATA_DIR/ALMANAC.CCLE_PC3.csv" \
	--synergy.train "CCLE.RPMI8226=$DATA_DIR/ALMANAC.CCLE_RPMI8226.csv" \
	--synergy.train "CCLE.SF268=$DATA_DIR/ALMANAC.CCLE_SF268.csv" \
	--synergy.train "CCLE.SF295=$DATA_DIR/ALMANAC.CCLE_SF295.csv" \
	--synergy.train "CCLE.SF539=$DATA_DIR/ALMANAC.CCLE_SF539.csv" \
	--synergy.train "CCLE.SKMEL28=$DATA_DIR/ALMANAC.CCLE_SKMEL28.csv" \
	--synergy.train "CCLE.SKMEL5=$DATA_DIR/ALMANAC.CCLE_SKMEL5.csv" \
	--synergy.train "CCLE.SKOV3=$DATA_DIR/ALMANAC.CCLE_SKOV3.csv" \
	--synergy.train "CCLE.SNB75=$DATA_DIR/ALMANAC.CCLE_SNB75.csv" \
	--synergy.train "CCLE.SR786=$DATA_DIR/ALMANAC.CCLE_SR786.csv" \
	--synergy.train "CCLE.SW620=$DATA_DIR/ALMANAC.CCLE_SW620.csv" \
	--synergy.train "CCLE.T47D=$DATA_DIR/ALMANAC.CCLE_T47D.csv" \
	--synergy.train "CCLE.U251MG=$DATA_DIR/ALMANAC.CCLE_U251MG.csv" \
	--synergy.train "CCLE.UACC257=$DATA_DIR/ALMANAC.CCLE_UACC257.csv" \
	--synergy.train "CCLE.UACC62=$DATA_DIR/ALMANAC.CCLE_UACC62.csv" \
	--synergy.train "CCLE.UO31=$DATA_DIR/ALMANAC.CCLE_UO31.csv" \
	--synergy.test "CCLE.A2058=$DATA_DIR/Merck.CCLE_A2058.csv" \
	--synergy.test "CCLE.A2780=$DATA_DIR/Merck.CCLE_A2780.csv" \
	--synergy.test "CCLE.A375=$DATA_DIR/Merck.CCLE_A375.csv" \
	--synergy.test "CCLE.A427=$DATA_DIR/Merck.CCLE_A427.csv" \
	--synergy.test "CCLE.CAOV3=$DATA_DIR/Merck.CCLE_CAOV3.csv" \
	--synergy.test "CCLE.COLO320=$DATA_DIR/Merck.CCLE_COLO320.csv" \
	--synergy.test "CCLE.DLD1=$DATA_DIR/Merck.CCLE_DLD1.csv" \
	--synergy.test "CCLE.EFM192A=$DATA_DIR/Merck.CCLE_EFM192A.csv" \
	--synergy.test "CCLE.ES2=$DATA_DIR/Merck.CCLE_ES2.csv" \
	--synergy.test "CCLE.HCT116=$DATA_DIR/Merck.CCLE_HCT116.csv" \
	--synergy.test "CCLE.HT144=$DATA_DIR/Merck.CCLE_HT144.csv" \
	--synergy.test "CCLE.HT29=$DATA_DIR/Merck.CCLE_HT29.csv" \
	--synergy.test "CCLE.KPL1=$DATA_DIR/Merck.CCLE_KPL1.csv" \
	--synergy.test "CCLE.LNCAPCLONEFGC=$DATA_DIR/Merck.CCLE_LNCAPCLONEFGC.csv" \
	--synergy.test "CCLE.LOVO=$DATA_DIR/Merck.CCLE_LOVO.csv" \
	--synergy.test "CCLE.MDAMB436=$DATA_DIR/Merck.CCLE_MDAMB436.csv" \
	--synergy.test "CCLE.MSTO211H=$DATA_DIR/Merck.CCLE_MSTO211H.csv" \
	--synergy.test "CCLE.NCIH1650=$DATA_DIR/Merck.CCLE_NCIH1650.csv" \
	--synergy.test "CCLE.NCIH2122=$DATA_DIR/Merck.CCLE_NCIH2122.csv" \
	--synergy.test "CCLE.NCIH23=$DATA_DIR/Merck.CCLE_NCIH23.csv" \
	--synergy.test "CCLE.NCIH460=$DATA_DIR/Merck.CCLE_NCIH460.csv" \
	--synergy.test "CCLE.NCIH520=$DATA_DIR/Merck.CCLE_NCIH520.csv" \
	--synergy.test "CCLE.NIHOVCAR3=$DATA_DIR/Merck.CCLE_NIHOVCAR3.csv" \
	--synergy.test "CCLE.PA1=$DATA_DIR/Merck.CCLE_PA1.csv" \
	--synergy.test "CCLE.OV90=$DATA_DIR/Merck.CCLE_OV90.csv" \
	--synergy.test "CCLE.RKO=$DATA_DIR/Merck.CCLE_RKO.csv" \
	--synergy.test "CCLE.RPMI7951=$DATA_DIR/Merck.CCLE_RPMI7951.csv" \
	--synergy.test "CCLE.SKMEL30=$DATA_DIR/Merck.CCLE_SKMEL30.csv" \
	--synergy.test "CCLE.SKMES1=$DATA_DIR/Merck.CCLE_SKMES1.csv" \
	--synergy.test "CCLE.SKOV3=$DATA_DIR/Merck.CCLE_SKOV3.csv" \
	--synergy.test "CCLE.SW620=$DATA_DIR/Merck.CCLE_SW620.csv" \
	--synergy.test "CCLE.SW837=$DATA_DIR/Merck.CCLE_SW837.csv" \
	--synergy.test "CCLE.T47D=$DATA_DIR/Merck.CCLE_T47D.csv" \
	--synergy.test "CCLE.UACC62=$DATA_DIR/Merck.CCLE_UACC62.csv" \
	--synergy.test "CCLE.VCAP=$DATA_DIR/Merck.CCLE_VCAP.csv" \
	--synergy.test "CCLE.ZR751=$DATA_DIR/Merck.CCLE_ZR751.csv"


