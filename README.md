## Investigation of Differential Gene Expression in Major Depressive Disorder Using Information Theory Supplemental information

### Code structure
* Main.py: Runs the main mutual information algorithm and stores results
    * mi.py: Contains the class for running mutual information algorithm
    * preprocessing.py: preprocesses the data into the correct format
    * process_reads.R: uses DeSeq2 to process raw read counts
* diffusion_map.py: generates diffusion_map plots
    * diffusion_map.R: Runs diffusion_map algorithm in destiny
    * plot_embeddings.py: plots the embeddings for various techniques (PCA, tsne, diffusionmap)
* gsea.py: runs gsea software
    * gsea_utils.py: Contains class for running gsea and building geneset database to run gsea
