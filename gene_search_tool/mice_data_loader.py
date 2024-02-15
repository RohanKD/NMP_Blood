if __name__ == "__main__":
    print("done")
    # figure out how to read .rds files that contain seurat objects
    # load data into an anndata dataframe and use same tool aas before to identify NMPs, just switch out data
    import rpy2
    #
    # result = pyreadr.read_r(
    #     'C:\\Users\\rohan\PycharmProjects\SciFaitEmory\code\gene_search_tool - Copy\seurat_object_E6.5.rds')  # also works for RData

    # done!
    # result is a dictionary where keys are the name of objects and the values python
    # objects. In the case of Rds there is only one object with None as key
    df = result[None]  # extract the pandas data frame
