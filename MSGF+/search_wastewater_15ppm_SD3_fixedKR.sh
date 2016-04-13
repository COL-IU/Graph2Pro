#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.hybrid.bp.fixedKR.mzid -d database/wastewater/hybrid_SD3-s2-63mer-k31-d1.contig-pep.fixedKR.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 2 -tda 0 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.hybrid.bp.fixedKR.mzid -showDecoy 1

#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.hybrid.frag.mzid -d database/wastewater/hybrid_SD3-s2-63mer-k31-d1.500.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 2 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.hybrid.frag.mzid -showDecoy 1

#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.mg.bp.fixedKR.mzid -d database/wastewater/MG_SRR1046369-SD3MG-s2-63mer-k31-d1.contig-pep.fixedKR.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 2 -tda 0 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.mg.bp.fixedKR.mzid -showDecoy 1

#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.mg.frag.mzid -d database/wastewater/MG_SRR1046369-SD3MG-s2-63mer-k31-d1.500.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 2 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.mg.frag.mzid -showDecoy 1

#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.tag.mzid -d database/wastewater/SD3-trimmed-k31-TAG.fgs.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.tag.mzid -showDecoy 1

#java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.db2.mzid -d database/wastewater/SD3_DBGraphPep2Pro_5.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
#java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.db2.mzid -showDecoy 1

java -Xmx32g -jar MSGFPlus.jar -s wastewater/SD3.mgf -o wastewater/SD3.hybrid.db2.mzid -d database/wastewater/SD3_hybrid_DBGraphPep2Pro_5.fasta -inst 1 -t 15ppm -ti -1,2 -mod Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1
java -Xmx16g -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i wastewater/SD3.hybrid.db2.mzid -showDecoy 1
