# login
ssh syang71@bayes.biostat.washington.edu

# upload files from local 
rsync -avz $[filepath] syang71@bayes.biostat.washington.edu:/home/users/syang71/Dataset/

# /home/users/syang71/Dataset/Larry_Dataset_normalized/stateFate_inVitro_normed_counts.mtx.gz"
# /home/users/syang71/SCLineage_ConstrativeLearning
# /home/users/syang71/kzlinlab/projects/scContrastiveLearn/git/scContrastiveLearn_Joshua


# copy 
cp scContrastiveLearning.py scContrastiveLearning_ver2.py

# submit job 
sbatch test-submit.slurm
squeue -u syang71

# check the size of a file 
ls -lh filename

# git push and etc 
git clone https://github.com/SZ-yang/SCSeq_LineageBarcoding
cd SCSeq_LineageBarcoding
git status    
git add Constrative_learning/*             or           git add .
git status
git commit -m 
git push origin master

git pull origin master





batch size 
10   1.08 220
15 1.47 
20 
# normalize 单个数据normlize ||x|| =1
# simclr 怎么normalize
# temperature
# test set 的 loss 