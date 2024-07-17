tdfread('Input/ImageInfo_List.txt');
fileID = fopen('Output/PearsonCoeff_Matlab_actual_scramble.txt','w');
fileID2 = fopen('Output/OverlapVolume_Matlab.txt','w');
fprintf(fileID,'FileName\tFrame\tAllele1\tAllele1_scramble\tAllele1_scramble_bootstrap\tAllele2\tAllele2_scramble\tAllele2_scramble_bootstrap\n');
fprintf(fileID2,'FileName\tAllele\tFrame\tNanog_Volume\tDNA_Volume\tRNA_Volume\tOverlap_Volume\n');
for i = 1:size(FileName,1)
    FileName(i,:)
    cnt=1;MatrixWithImages=[];
    for t = 1:NumberFrames(i)
        for z = 1:NumberSlices(i)
            for c = 1:NumberChannels(i)
                MatrixWithImages(:,:,c,z,t) = imread(strtrim(FileName(i,:)),cnt);
                cnt=cnt+1;
            end
        end
    end
    for t = 1:NumberFrames(i)
        nanog1=[];dna1=[];nanog2=[];dna2=[];ovlp_vol1=0;nanog_vol1=0;
        dna_vol1=0;rna_vol1=0;ovlp_vol2=0;nanog_vol2=0;dna_vol2=0;rna_vol2=0;
        for z = 1:NumberSlices(i)
            nanogz = double(MatrixWithImages(:,:,1,z,t));
            dnaz = double(MatrixWithImages(:,:,2,z,t));
            
            if(NumberAlleles(i)>1)
                nanogz1_mask = double(MatrixWithImages(:,:,4,z,t)>0);
                dnaz1_mask = double(MatrixWithImages(:,:,6,z,t)>0);
                rnaz1_mask = double(MatrixWithImages(:,:,8,z,t)>0);
                indz1 = find(dnaz1_mask>0);
                nanogz1 = nanogz(indz1);
                dnaz1 = dnaz(indz1);
                nanog1=[nanog1;nanogz1];dna1=[dna1;dnaz1];
                ovlp_vol1=ovlp_vol1+size(find(dnaz1_mask>0 & nanogz1_mask>0),1);
                nanog_vol1 = nanog_vol1 + size(find(nanogz1_mask>0),1);
                dna_vol1= dna_vol1+ size(find(dnaz1_mask>0),1);
                rna_vol1= rna_vol1+ size(find(rnaz1_mask>0),1);

                nanogz2_mask = double(MatrixWithImages(:,:,5,z,t)>0);
                dnaz2_mask = double(MatrixWithImages(:,:,7,z,t)>0);
                rnaz2_mask = double(MatrixWithImages(:,:,9,z,t)>0);
                indz2 = find(dnaz2_mask>0);
                nanogz2 = nanogz(indz2);
                dnaz2 = dnaz(indz2);
                nanog2=[nanog2;nanogz2];dna2=[dna2;dnaz2];
                ovlp_vol2=ovlp_vol2+size(find(dnaz2_mask>0 & nanogz2_mask>0),1);
                nanog_vol2 = nanog_vol2 + size(find(nanogz2_mask>0),1);
                dna_vol2= dna_vol2+ size(find(dnaz2_mask>0),1);
                rna_vol2= rna_vol2+ size(find(rnaz2_mask>0),1);
            else
                nanogz1_mask = double(MatrixWithImages(:,:,4,z,t)>0);
                dnaz1_mask = double(MatrixWithImages(:,:,5,z,t)>0);
                rnaz1_mask = double(MatrixWithImages(:,:,6,z,t)>0);
                indz1 = find(dnaz1_mask>0);
                nanogz1 = nanogz(indz1);
                dnaz1 = dnaz(indz1);
                nanog1=[nanog1;nanogz1];dna1=[dna1;dnaz1];
                ovlp_vol1=ovlp_vol1+size(find(dnaz1_mask>0 & nanogz1_mask>0),1);
                nanog_vol1 = nanog_vol1 + size(find(nanogz1_mask>0),1);
                dna_vol1= dna_vol1+ size(find(dnaz1_mask>0),1);
                rna_vol1= rna_vol1+ size(find(rnaz1_mask>0),1);
            end
        end
        if(size(nanog1,1)>1)
            [p1,r1]=corrcoef(nanog1,dna1);p1_rnd=[];
            for bt=1:10
                nanog1_rnd = nanog1(randperm(size(nanog1,1)),:);
                [p1_rnd_bt,~]=corrcoef(nanog1_rnd,dna1);
                p1_rnd=[p1_rnd;p1_rnd_bt(1,2)];
            end
            pcorr1=p1(1,2);pcorr_rnd1= p1_rnd_bt(1,2); pcorr_mean_rnd1=mean(p1_rnd);
        else
            pcorr1=NaN;pcorr_rnd1=NaN;pcorr_mean_rnd1=NaN;
        end
        if(NumberAlleles(i)>1 && size(nanog2,1)>1)
            [p2,r2]=corrcoef(nanog2,dna2);p2_rnd=[];            
            for bt=1:10
                nanog2_rnd = nanog2(randperm(size(nanog2,1)),:);
                [p2_rnd_bt,~]=corrcoef(nanog2_rnd,dna2);
                p2_rnd=[p2_rnd;p2_rnd_bt(1,2)];
            end
            pcorr2=p2(1,2);pcorr_rnd2= p2_rnd_bt(1,2); pcorr_mean_rnd2=mean(p2_rnd);
        else
            pcorr2=NaN;pcorr_rnd2=NaN;pcorr_mean_rnd2=NaN;
        end
        fprintf(fileID2,'%s\t1\t%d\t%d\t%d\t%d\t%d\n',FileName(i,:),t,nanog_vol1,dna_vol1,rna_vol1,ovlp_vol1);
        if(NumberAlleles(i)>1)
            fprintf(fileID2,'%s\t2\t%d\t%d\t%d\t%d\t%d\n',FileName(i,:),t,nanog_vol2,dna_vol2,rna_vol2,ovlp_vol2);
            fprintf(fileID,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',FileName(i,:),t,pcorr1,pcorr_rnd1,pcorr_mean_rnd1,pcorr2,pcorr_rnd2,pcorr_mean_rnd2);
        else
            fprintf(fileID,'%s\t%d\t%d\t%d\t%d\n',FileName(i,:),t,pcorr1,pcorr_rnd1,pcorr_mean_rnd1);
        end       
    end
    sprintf('%s',FileName(i,:))
end