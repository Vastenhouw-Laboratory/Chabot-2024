#!/usr/bin/perl
open F0, "Input/Filelist_ShortestDist_MCP.txt" or die;
#12.01_1k_9_1_Nanog_Shortest_Distance_to_Surfaces_Surfaces=MCP.csv
while($fname=<F0>){
	chomp $fname;
	@tmp = split "_Nanog",$fname;$nuclei_id=$tmp[0];
	push @nuclei_list,$nuclei_id;
	$fname1 = $nuclei_id."_Nanog_Diameter.csv";
	open F1, "Input/Statistics/$fname1" or die;
	#Diameter X,Diameter Y,Diameter Z,Unit,Category,Collection,Time,TrackID,ID,
	#0.688011,0.688011,1.37602,µm,Spot,Diameter,1,1000000000,0,
	$line=<F1>;$line=<F1>;$line=<F1>;$line=<F1>;%cutoff=();
	while($line=<F1>){
		chomp $line;
		@spl = split ",", $line;
		$cutoff{$spl[$#spl-1]}=$spl[0];
	}
	open F1, "Input/Statistics/$fname" or die;
	#Shortest Distance to Surfaces,Unit,Category,Surfaces,Time,TrackID,ID,
	#-0.286843,,Spot,MCP,1,1000000000,0,
	$line=<F1>;$line=<F1>;$line=<F1>;$line=<F1>;
	while($line=<F1>){
		chomp $line;
		@spl = split ",", $line;
		$trackid=$spl[$#spl-2];$time=$spl[$#spl-3];
		$spotnum{$nuclei_id}{$trackid}{$time}++;
		if($spl[0]<$cutoff{$spl[$#spl-1]}){$ovlp_MCP{$nuclei_id}{$trackid}{$time}++;}
	}
}
open F1, ">Output/AllNuclei_TrackID_FusionTimes_NumNanogClust-Ovlp-Cutoff_RadiusZ.txt" or die;
print F1 "\t\t\tNanogClustOvlpFraction_TimeFrames\n";
print F1 "Nuclei_ID\tTrackID\tFusionTimes\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\n";
for($i=0;$i<scalar @nuclei_list;$i++){
	$nuclei_id = $nuclei_list[$i];
	foreach $trckId(keys %{$spotnum{$nuclei_id}}){
		print F1 "$nuclei_id\t$trckId\t";
		$line="";
		$time_fusion{$nuclei_id}{$trckId}=-1;
		$prv_spotnum=-1;
		for($t=1;$t<=20;$t++){
			if(! exists $spotnum{$nuclei_id}{$trckId}{$t}){
				$spotnum{$nuclei_id}{$trckId}{$t}=0;
			}
			if(! exists $ovlp_MCP{$nuclei_id}{$trckId}{$t}){
				$ovlp_MCP{$nuclei_id}{$trckId}{$t}=0;
			}
			$line.="$ovlp_MCP{$nuclei_id}{$trckId}{$t}/$spotnum{$nuclei_id}{$trckId}{$t}\t";
			if($prv_spotnum >1 && $spotnum{$nuclei_id}{$trckId}{$t} ==1){
				if($time_fusion{$nuclei_id}{$trckId}==-1){$time_fusion{$nuclei_id}{$trckId}=$t;}
				else{$time_fusion{$nuclei_id}{$trckId}.="|".$t;}
			}
			if($spotnum{$nuclei_id}{$trckId}{$t}>0){$prv_spotnum=$spotnum{$nuclei_id}{$trckId}{$t};}
		}
		print F1 "$time_fusion{$nuclei_id}{$trckId}\t$line\n";
	}
}
