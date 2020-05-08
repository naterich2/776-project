def write_gct(dataframe,outfile):
    """Write a GCT file from a dataframe

    Parameters
    ----------
    dataframe : sampleXcolumn dataframe
    outfile : file to write to

    Returns
    -------
    TODO

    """
    dataframe = dataframe.T
    dataframe.index.name = 'NAME'
    dataframe.insert(0,column='Description',value='')
    with open(outfile,'w') as out:
        print('#1.2',file=out)
        print(dataframe.shape[0],(dataframe.shape[1]-1),file=out)
    dataframe.to_csv(outfile,sep='\t',mode='a',header=True)

def write_cls(meta,outfile):
    """TODO: Docstring for write_cls.

    Parameters
    ----------
    meta : TODO
    outfile : TODO

    Returns
    -------
    TODO

    """
    pass
