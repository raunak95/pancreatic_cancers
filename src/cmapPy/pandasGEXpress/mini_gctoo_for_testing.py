"""
Creates a small GCToo instance (with representative examples of typically found fields); can use for testing.

ex:
    import mini_gctoo_for testing
    my_mini_gctoo = mini_gctoo_for_testing.make()
"""
import logging
import pandas
import numpy
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger

__author__ = 'Oana Enache'
__email__ = 'oana@broadinstitute.org'

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def make(convert_neg_666=True):
    """
    Creates a small GCToo instance (with representative examples of typically found fields); can use for testing.
    """
    # metadata examples; should be one of each type reasonable to find
    id_vals = ["LJP007_MCF10A_24H:TRT_CP:BRD-K93918653:3.33", "MISC003_A375_24H:TRT_CP:BRD-K93918653:3.33",
               "LJP007_MCF7_24H:TRT_POSCON:BRD-K81418486:10", "LJP007_MCF7_24H:TRT_POSCON:BRD-A61304759:10",
               "LJP007_MCF7_24H:CTL_VEHICLE:DMSO:-666", "LJP007_MCF7_24H:TRT_CP:BRD-K64857848:10"]
    count_cv = ["14|15|14", "13|14|13",
                "13|15|14|14|15|14|14|13|14|15|15|14|14|15|14|15|14|14|15|14|15|14|14|14|14|14|14|15|14|14|15|14|14|14|14|13|14|14|14|14|14|14|15|14|13|13|15|14|14|15|14|14|14|15|13|13|15|13|14|13|13|14|14|14|14|13",
                "13", "13", "14"]
    distil_ss = [9.822065353, 6.8915205, 1.35840559, 5.548898697, 3.355231762, 4.837643147]
    zmad_ref = ["population", "population", "population", "population", "population", "population"]
    distil_nsample = [3, 3, 66, 2, 9, 111111]
    mfc_plate_id = ["-666", "-666", "-666", "-666", "-666", "-666"]

    # build metadata dataframe
    mini_meta_dict = {}
    mini_meta_dict["id"] = id_vals
    mini_meta_dict["count_cv"] = count_cv
    mini_meta_dict["distil_ss"] = distil_ss
    mini_meta_dict["zmad_ref"] = zmad_ref
    mini_meta_dict["distil_nsample"] = distil_nsample
    mini_meta_dict["mfc_plate_id"] = mfc_plate_id
    mini_row_metadata = pandas.DataFrame(mini_meta_dict,
                                         columns=['id', 'count_cv', 'distil_nsample', 'distil_ss', 'mfc_plate_id', 'zmad_ref'])

    if convert_neg_666:
        mini_row_metadata = mini_row_metadata.replace([-666, "-666", -666.0], [numpy.nan, numpy.nan, numpy.nan])
        # if all values in a column are nanpandas.Series(mini_row_metadata.isna().sum() == mini_row_metadata.shape[0]) convert dtype of that column to float
        all_nan_columns = (mini_row_metadata.isnull().sum() == numpy.array(mini_row_metadata.shape[0])).nonzero()[0]
        mini_row_metadata = mini_row_metadata.astype({d: 'float' for d in mini_row_metadata.columns[all_nan_columns.tolist()]})
    else:
        mini_row_metadata = mini_row_metadata.replace([-666, -666.0], ["-666", "-666"])

    # for now (at least) col and row metadata are the same
    mini_col_metadata = mini_row_metadata.copy()

    # data example values
    r1 = [1, 2, 3, 4, 5, 6]
    r2 = [4.3, 4.5, 4.3, 4.3, 4.3, 4.3]
    r3 = [7, 8, 9, 0, 1.23476, 9.758320]
    r4 = [0.11, 3.3456356, 2.345667, 9.822065353, 4.78865099, 4.7886]
    r5 = [-0.11, -3.3456356, -2.345667, -9.822065353, -4.78865099, -4.7886]
    r6 = [1, -2, 3, -4, 5, -6]

    # build data dataframe
    mini_data_mat = pandas.DataFrame([r1, r2, r3, r4, r5, r6], dtype=numpy.float32)
    mini_data_mat.index = id_vals
    mini_data_mat.columns = id_vals

    # instantiate & assign attributes of GCToo instance
    mini_version = "GCTX1.0"
    mini_src = "mini_gctoo.gctx"

    mini_row_metadata_df = mini_row_metadata
    mini_row_metadata_df.set_index("id", inplace=True, drop=True)
    mini_row_metadata.index.name = "rid"
    mini_row_metadata_df.columns.name = "rhd"

    mini_col_metadata_df = mini_col_metadata
    mini_col_metadata_df.set_index("id", inplace=True, drop=True)
    mini_col_metadata.index.name = "cid"
    mini_col_metadata_df.columns.name = "chd"

    mini_data_df = mini_data_mat
    mini_data_df.index.name = "rid"
    mini_data_df.columns.name = "cid"

    logger.debug("Making mini_gctoo instance...")
    mini_gctoo = GCToo.GCToo(data_df=mini_data_df, row_metadata_df=mini_row_metadata_df,
                             col_metadata_df=mini_col_metadata_df, src=mini_src, version=mini_version)

    return mini_gctoo
