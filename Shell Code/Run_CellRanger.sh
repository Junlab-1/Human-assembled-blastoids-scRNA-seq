# !/bin/bash
# Run Blastoids, ICM, and TE spheroid samples(human samples). 
nohup cellranger count --id=Blastoid --transcriptome=/home/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW585/ --sample=Blastoid --create-bam=false --include-introns=true > Blastoids.txt 2>&1 &
nohup cellranger count --id=ICM --transcriptome=/home/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW587/ --sample=ICM --create-bam=false --include-introns=true > ICM.txt 2>&1 &
nohup cellranger count --id=TE_spheroid --transcriptome=/home/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW586/ --sample=TE_spheroid --create-bam=false --include-introns=true > TE_spheroid.txt 2>&1 &
# Run Chimeric Blastoids and Mouse ICM data
nohup cellranger count --localcores=16 --id=Chimeric_Blastoids_mouse --transcriptome=/home/reference/cellranger_ref/mouse/refdata-gex-mm10-2020-A --fastqs=/home/data/wuhao/scdata/chimeric_scRNA-seq/01.RawData/LW817 --sample=LW817-SI_TT_A10_22V5CYLT4,LW817-SI_TT_A10_22V5YHLT4,LW817-SI_TT_A10_22V5YJLT4 --create-bam=true --include-introns=true > ChimericB_mouse.txt 2>&1 &
nohup cellranger count --localcores=16 --id=Chimeric_Blastoids_human --transcriptome=/home/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=/home/data/wuhao/scdata/chimeric_scRNA-seq/01.RawData/LW817 --sample=LW817-SI_TT_A10_22V5CYLT4,LW817-SI_TT_A10_22V5YHLT4,LW817-SI_TT_A10_22V5YJLT4 --create-bam=true --include-introns=true > ChimericB_human.txt 2>&1 &
nohup cellranger count --localcores=16 --id=MouseICM --transcriptome=/home/reference/cellranger_ref/mouse/refdata-gex-mm10-2020-A --fastqs=/home/DX6/jwulab/S227384/data/wuhao/scdata/chimeric_scRNA-seq/01.RawData/LW818 --sample=LW818-SI_TT_A11_22V5CYLT4,LW818-SI_TT_A11_22V5YHLT4,LW818-SI_TT_A11_22V5YJLT4 --create-bam=true --include-introns=true > MouseICM.txt 2>&1 &

