#!/bin/sh

# When running this script via a scheduler, you will need to specify the starting directory.
# To do this, you can add an explict "cd <starting directory>" here. This starting directory is
# assumed to be the directory that contains the "predict_synergy" program.

DRUG_FEATURES=../preprocess/ALMANAC_and_Merck.obable.FP2.csv
CELL_FEATURES=../preprocess/CCLE_DepMap_18q3_RNAseq_RPKM.LINCS1000.tsv

# Replace "output" with the name of the directory that you have created to store the output files
OUTPUT_DIR=output

DATA_DIR=../preprocess/data_bliss
#DATA_DIR=../preprocess/data_loewe

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

cv_predict()
{

	# Argument $1 is the output filename
	# Optional argument two is either the drug or cell features
	# Optional argument three is either the drug or cell features
	
	mpirun -np $NUM_NODES \
		--oversubscribe \
		--bind-to socket \
		-x OMP_NUM_THREADS \
		./predic_synergy \
		--enable-pair-prediction \
		--forest.size $FOREST_SIZE \
		-o $1$2$3 \
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
		--synergy.train "CCLE.UO31=$DATA_DIR/ALMANAC.CCLE_UO31.csv"
}

# Add a space before the optional second and third function arguments
for ITER in {1..25}
do
	cv_predict "$OUTPUT_DIR/cell.drug.$ITER.$FOREST_SIZE.out" " --drug $DRUG_FEATURES" " --cell $CELL_FEATURES"
	#cv_predict "$OUTPUT_DIR/cell.hot-drug.$ITER.$FOREST_SIZE.out" " --hot.drug" " --cell $CELL_FEATURES"
	cv_predict "$OUTPUT_DIR/cell.no-drug.$ITER.$FOREST_SIZE.out" " --cell $CELL_FEATURES"
	
	cv_predict "$OUTPUT_DIR/hot-tissue.drug.$ITER.$FOREST_SIZE.out" " --drug $DRUG_FEATURES" " --hot.tissue"
	#cv_predict "$OUTPUT_DIR/hot-tissue.hot-drug.$ITER.$FOREST_SIZE.out" " --hot.drug" " --hot.tissue"
	cv_predict "$OUTPUT_DIR/hot-tissue.no-drug.$ITER.$FOREST_SIZE.out" " --hot.tissue"
	
	#cv_predict "$OUTPUT_DIR/hot-cell.drug.$ITER.$FOREST_SIZE.out" " --drug $DRUG_FEATURES" " --hot.cell"
	#cv_predict "$OUTPUT_DIR/hot-cell.hot-drug.$ITER.$FOREST_SIZE.out" " --hot.drug" " --hot.cell"
	#cv_predict "$OUTPUT_DIR/hot-cell.no-drug.$ITER.$FOREST_SIZE.out" " --hot.cell"
	
	cv_predict "$OUTPUT_DIR/no-cell.drug.$ITER.$FOREST_SIZE.out" " --drug $DRUG_FEATURES"
	#cv_predict "$OUTPUT_DIR/no-cell.hot-drug.$ITER.$FOREST_SIZE.out" " --hot.drug"
	cv_predict "$OUTPUT_DIR/no-cell.no-drug.$ITER.$FOREST_SIZE.out"
done

