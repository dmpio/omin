import omin

def compare(experiment_object=None, fraction=None, genotype=None, treatment=None,comparison= None):
    """Returns a tuple containing the Log2 fold change and p-values for a given comparison of an Experiment object.

    Parameters
    ----------
    experiment_object : (:obj)
    fraction : str
    genotype : str or list
        Use a string when comparison = 'within'. Use a two element list when comparison = 'between'.
    treatment : str or list
        Use a string when comparison = 'between'. When comparison = 'with' use a two element list with the
        first element is the perturbation and the second is the control.
    comparison : str
        Can either be 'within' for a within genotype comparison or 'between' for a between genotype comparison.

    Returns
    -------
    (lfc,pval) : tuple(DataFrame,DataFrame)
        First element is a DataFrame of the Log2 fold changes, second element is a DataFrame of the p-values.

    Examples
    --------
    Making a within genotype comparison;
    >>>lfc,pval = compare(experiment_object,fraction="Acetyl",genotype=["ko","wt"],treatment="nonex",comparison="between")
    >>>omin.vis.plotByMito(lfc,pval,mitodex,nonmitodex)

    Making a between genotype comparison;
    >>>lfc,pval = compare(experiment_object,fraction="Phospho",genotype="ko",treatment=["immediate_post","nonex"],comparison="within")
    >>>omin.vis.plotByMito(lfc,pval,mitodex,nonmitodex)


    """
    if not isinstance(experiment_object,omin.Experiment):
        raise NotImplementedError("raw_file must be a omin.Experiment object.")

    if comparison is None:
        raise NotImplementedError("Please enter one of the these: 'within' or 'between'")

    if comparison.lower() == "within":
        if type(treatment) is not list:
            raise NotImplementedError("The parameter 'treatment' must be a list for a within genotype comparison.")
        if not isinstance(genotype,str):
            raise NotImplementedError("The parameter 'genotype' must be a string for a within genotype comparison")
        #Make the comparison
        within_geno = tuple(experiment_object.peptides.__dict__[fraction].__dict__[genotype].__dict__[treat] for treat in treatment)
        # Create a descriptive column name
        col_name = " ".join([fraction,genotype,treatment[0]+"/"+treatment[1]])
        #ttest
        pval = omin.ttester(within_geno[0],within_geno[1],new_column_name = col_name)
        #Log2 fold change
        lfc = omin.log2FC(within_geno[0],within_geno[1],new_column_name = col_name)
        return lfc,pval

    if comparison.lower() == "between":
        if type(genotype) is not list:
            raise NotImplementedError("The parameter 'genotype' must be a list of two genotypes for between genotype comparison.")
        if type(treatment) is not str:
            raise NotImplementedError("The paramter 'treatment' must be a string for between genotype comparison.")
        between_geno = tuple(experiment_object.peptides.__dict__[fraction].__dict__[geno].__dict__[treatment] for geno in genotype)

        # Create a descriptive column name
        col_name = " ".join([fraction,"/".join(genotype)])
        pval = omin.ttester(between_geno[0],between_geno[1],new_column_name = col_name)
        lfc = omin.log2FC(between_geno[0],between_geno[1],new_column_name = col_name)
        return lfc,pval
