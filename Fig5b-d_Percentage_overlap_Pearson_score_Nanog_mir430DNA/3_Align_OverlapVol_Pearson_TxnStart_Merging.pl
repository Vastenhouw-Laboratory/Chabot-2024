#!/usr/bin/perl
open F1, "Output/OverlapVolume_Matlab.txt" or die;
#FileName        Allele  Frame   Nanog_Volume    DNA_Volume      RNA_Volume      Overlap_Volume
#Input/Nanog_miR430_Nuclei_Tif_Format/09.01_1k_3_1.ims.tif     1       1       0       57      0       0
$line=<F1>;
while($line=<F1>){
	chomp $line;
	@spl = split "\t", $line;
	@tmp = split "/",$spl[0];
	$tmp[$#tmp]=~ s/\.ims\.tif//;$tmp[$#tmp]=~s/ //g;
	$ID=$tmp[$#tmp]."_".$spl[1];
	$Nanog_Vol{$ID}{$spl[2]}=$spl[3];
	$DNA_Vol{$ID}{$spl[2]}=$spl[4];
	$RNA_Vol{$ID}{$spl[2]}=$spl[5];
	if($spl[3]>0 and $spl[4]>0){$Ovlp_Vol_Nanog_DNA{$ID}{$spl[2]}=$spl[6];}
	else{$Ovlp_Vol_Nanog_DNA{$ID}{$spl[2]}="NA";}
}
open F1, "Output/PearsonCoeff_Matlab_actual_scramble.txt" or die;
#FileName        Frame   Allele1 Allele1_scramble        Allele1_scramble_bootstrap      Allele2 Allele2_scramble        Allele2_scramble_bootstrap
#Input/Nanog_miR430_Nuclei_Tif_Format/09.01_1k_3_1.ims.tif     1       -1.876896e-01   2.848990e-01    2.391612e-02    2.656201e-01    1.172360e-01    2.681598e-02
$line=<F1>;
while($line=<F1>){
	chomp $line;
	@spl = split "\t", $line;
	@tmp = split "/",$spl[0];
	$tmp[$#tmp]=~ s/\.ims\.tif//;$tmp[$#tmp]=~s/ //g;
	$ID1=$tmp[$#tmp]."_1";
	$Pearson_act{$ID1}{$spl[1]}=$spl[2];
	$Pearson_scr{$ID1}{$spl[1]}=$spl[3];
	if($spl[5] ne ""){
		$ID2=$tmp[$#tmp]."_2";
		$Pearson_act{$ID2}{$spl[1]}=$spl[5];
		$Pearson_scr{$ID2}{$spl[1]}=$spl[6];
	}
}
open F1, "Output/AllNuclei_TrackID_FusionTimes_NumNanogClust-Ovlp-Cutoff_RadiusZ.txt" or die;
#			NanogClustOvlpFraction_TimeFrames
#Nuclei_ID	TrackID	FusionTimes	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
#09.01_1k_3_1	1000000003	11	0/0	0/0	0/0	1/1	1/1	2/2	2/2	2/2	2/2	2/2	1/1	2/2	2/2	2/2	2/2	2/2	3/3	2/2	2/2	2/2
$line=<F1>;$line=<F1>;
while($line=<F1>){
	chomp $line;
	@spl = split "\t", $line;
	$ID=$spl[0]."_".$spl[1];
	if($spl[2] ne -1){$FusionTimes{$ID}=$spl[2];}
	$OvlpClustFrc{$ID}=join("\t",@spl[3..$#spl]);
}

open F1, "Input/MasterTable_AllNuclei_AlleleNum_TrackID_TxnTime_30Oct2023.csv" or die;
#Nuclei_ID,AlleleNum,TrackID,Time_Activation
#09.01_1k_3_1,2,1000000003,9
open F0, ">Output/AllNuclei_Allele_All_Values_Pooled.txt" or die;
print F0 "ID\tTrackID\tMerge_Flag\tTime\tTxn_Time\tMerge_Time\tRelTime_Txn\tRelTime_Fus\tNumNanogClust\tNumNanogClustOvlpDNA\tNanog_Volume\tDNA_Volume\tRNA_Volume\tOverlap_Volume\tPearsonCoeff\tPearsonCoeff_scramble\n";
$line=<F1>;
while($line=<F1>){
	chomp $line;
	@spl = split ",", $line;
	$ID=$spl[0]."_".$spl[1];$ID1=$spl[0]."_".$spl[2];
	$txn_time=$spl[3];
	@fus_time=split '\|',$FusionTimes{$ID1};
	@clustovlp=split "\t",$OvlpClustFrc{$ID1};
	$more_thn_1_clust_bef_txn=0;
	for($i=0;$i<scalar @clustovlp; $i++){
		@tmp = split "/", $clustovlp[$i];
		if($tmp[1] >1 and (($i+1) <$txn_time)){$more_thn_1_clust_bef_txn=1;} 
	}
	$smallest_pos=100;$largest_neg=-100;
	for($i=0;$i< scalar @fus_time;$i++){
		$d=$txn_time-$fus_time[$i];
		if($d>=0){if($d < $smallest_pos){$smallest_pos=$d;}}
		if($d<0){if($largest_neg < $d){$largest_neg=$d;}}
	}
	if($more_thn_1_clust_bef_txn > 0){
		if($smallest_pos != 100){
			$time_merge = $txn_time - $smallest_pos;
		}else{
			if($largest_neg >= -3){$time_merge = $txn_time - $largest_neg;}
			else{$time_merge="NA";}
		}
	}else{$time_merge="NA";}
	if($time_merge ne "NA"){$flag_merge=1;}else{$flag_merge=0;}
	$time = 0;
        for($i=0;$i<scalar @clustovlp;$i++){
                $time++;$time_rel_txn=$time-$txn_time;
                if($flag_merge !=0){$time_rel_merge=$time-$time_merge;}
                else{$time_rel_merge="NA";}
                @tmp = split "/",$clustovlp[$i];
                $NumNanogClust=$tmp[1];
                if($Nanog_Vol{$ID}{$time}>0 and $DNA_Vol{$ID}{$time}>0){$NumNanogClustOvlp=$tmp[0];}
                else{$NumNanogClustOvlp="NA";}
                print F0 "$ID\t$spl[2]\t$flag_merge\t$time\t$txn_time\t$time_merge\t$time_rel_txn\t$time_rel_merge\t$NumNanogClust\t$NumNanogClustOvlp\t$Nanog_Vol{$ID}{$time}\t$DNA_Vol{$ID}{$time}\t$RNA_Vol{$ID}{$time}\t$Ovlp_Vol_Nanog_DNA{$ID}{$time}\t$Pearson_act{$ID}{$time}\t$Pearson_scr{$ID}{$time}\n";
        }	
}	

