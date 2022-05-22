# documentation
# pydamage.readthedocs.io

#git link
# https://github.com/maxibor/pydamage


#citation
#site the bioarchive and the git in the same way as they did in thier paper 

#start install perceedure
#if in base when ~/src than diactivate and follow with the clone and installation.
conda deactivate


######

cd ~/src


#copy the url for the git clone in the future. 
git clone https://github.com/maxibor/pydamage
cd pydamage
git checkout dev
conda env create -f environment.yml
conda activate pydamage
pip install -e .
#then test it
pydamage --help

#run 
# pydamage analyze aligned.bam
#file is in data
pydamage analyze ERR3678595.aligned.sorted.bam

pydamage.main.pydamage_analyze_group(ERR3678595.aligned.sorted.bam, wlen=30, show_al=False, process=1, outdir='', plot=False, verbose=False, force=False)

# step 1, recreate from the raw data, or used what they dcreated a bam file, and than run the  pydamage in herelike in the git

# fastq dump the raw data in create bam file and than apply to the proceedure. 