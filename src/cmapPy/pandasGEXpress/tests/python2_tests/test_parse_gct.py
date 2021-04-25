import unittest
import logging
import os
import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger
import cmapPy.pandasGEXpress.parse_gct as pg
import cmapPy.pandasGEXpress.GCToo as GCToo


FUNCTIONAL_TESTS_PATH = "cmapPy/pandasGEXpress/tests/functional_tests/"

logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestParseGct(unittest.TestCase):
    def test_read_version_and_dims(self):
        ### v1.3 case
        version1 = "1.3"
        version1_as_string = "GCT1.3"
        dims1 = ["10", "15", "3", "4"]
        fname1 = "testing_testing1"

        f1 = open(fname1, "w")
        f1.write(("#" + version1 + "\n"))
        f1.write((dims1[0] + "\t" + dims1[1] + "\t" + dims1[2] + "\t" + dims1[3] + "\n"))
        f1.close()

        (actual_version, n_rows, n_cols, n_rhd, n_chd) = pg.read_version_and_dims(fname1)
        self.assertEqual(actual_version, version1_as_string)
        self.assertEqual(n_rows, int(dims1[0]))
        self.assertEqual(n_chd, int(dims1[3]))

        # Remove the file created
        os.remove(fname1)

        ### v1.2 case
        version2 = "1.2"
        version2_as_string = "GCT1.2"
        dims2 = ["10", "15"]
        fname2 = "testing_testing2"

        f2 = open(fname2, "w")
        f2.write(("#" + version2 + "\n"))
        f2.write((dims2[0] + "\t" + dims2[1] + "\n"))
        f2.close()

        (actual_version, n_rows, n_cols, n_rhd, n_chd) = pg.read_version_and_dims(fname2)
        self.assertEqual(actual_version, version2_as_string)
        self.assertEqual(n_rows, int(dims2[0]))
        self.assertEqual(n_cols, int(dims2[1]))

        # Remove the file created
        os.remove(fname2)

    def test_parse_into_3_df(self):
        gct_filepath = os.path.join(FUNCTIONAL_TESTS_PATH, "test_l1000.gct")
        e_dims = [978, 377, 11, 35]
        (row_metadata, col_metadata, data) = pg.parse_into_3_df(
            gct_filepath, e_dims[0], e_dims[1], e_dims[2], e_dims[3], None)

        # Check shapes of outputs
        self.assertTrue(row_metadata.shape == (e_dims[0], e_dims[2]),
                        ("row_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(row_metadata.shape,
                                                                  e_dims[0], e_dims[2]))
        self.assertTrue(col_metadata.shape == (e_dims[1], e_dims[3]),
                        ("col_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(col_metadata.shape,
                                                                  e_dims[1], e_dims[3]))
        self.assertTrue(data.shape == (e_dims[0], e_dims[1]),
                        ("data.shape = {} " +
                         "but expected it to be ({}, {})").format(data.shape,
                                                                  e_dims[0], e_dims[1]))

        # Type-check the data
        self.assertTrue(isinstance(data.iloc[0, 0], np.float32), "The data should be a float32, not {}".format(type(data.iloc[0, 0])))

        # Check a few values
        correct_val = np.float32(11.3818998337)
        self.assertTrue(data.iloc[0, 0] == correct_val,
                        ("The first value in the data matrix should be " +
                         "{} not {}").format(correct_val, data.iloc[0, 0]))
        correct_val = np.float32(5.1256)
        self.assertTrue(data.iloc[e_dims[0] - 1, e_dims[1] - 1] == correct_val,
                        ("The last value in the data matrix should be " +
                         str(correct_val) + " not {}").format(data.iloc[e_dims[0] - 1, e_dims[1] - 1]))
        correct_str = "LUA-4000"
        self.assertTrue(row_metadata.iloc[2, 3] == correct_str,
                        ("The 3rd row, 4th column of the row metadata should be " +
                         correct_str + " not {}").format(row_metadata.iloc[2, 3]))
        correct_str = 57
        self.assertTrue(col_metadata.iloc[e_dims[1] - 1, 0] == correct_str,
                        ("The last value in the first column of column metadata should be " +
                         str(correct_str) + " not {}").format(col_metadata.iloc[e_dims[1] - 1, 0]))

        # Check headers
        correct_str = "LJP005_A375_24H_X1_B19:P24"
        self.assertTrue(col_metadata.index.values[e_dims[1] - 1] == correct_str,
                        ("The last column metadata index should be " +
                         correct_str + " not {}").format(col_metadata.index.values[e_dims[1] - 1]))
        correct_str = "bead_batch"
        self.assertTrue(list(col_metadata)[3] == correct_str,
                        ("The fourth column metadata index value should be " +
                         correct_str + " not {}").format(list(col_metadata)[3]))
        correct_str = "203897_at"
        self.assertTrue(row_metadata.index.values[e_dims[0] - 1] == correct_str,
                        ("The last row metadata index value should be " + correct_str +
                         " not {}").format(row_metadata.index.values[e_dims[0] - 1]))
        self.assertTrue(data.index.values[e_dims[0] - 1] == correct_str,
                        ("The last data index value should be " + correct_str +
                         " not {}").format(data.index.values[e_dims[0] - 1]))

    def test_assemble_row_metadata(self):
        #simple happy path
        full_df = pd.DataFrame(
            [["id", "rhd1", "id", "cid1", "cid2"],
             ["chd1", "", "", "a", "b"],
             ["chd2", "", "", "55", "61"],
             ["chd3", "", "", "nah", "nope"],
             ["rid1", "C", "1.0", "0.3", "0.2"],
             ["rid2", "D", "2.0", np.nan, "0.9"]])
        full_df_dims = [2, 2, 2, 3]
        row_index = pd.Index(["rid1", "rid2"], name="rid")
        col_index = pd.Index(["rhd1", "id"], name="rhd")
        e_row_df = pd.DataFrame([["C", 1.0], ["D", 2.0]],
                                index = row_index,
                                columns = col_index)
        row_df = pg.assemble_row_metadata(full_df, full_df_dims[3],
                                          full_df_dims[0], full_df_dims[2])
        logger.debug("simple happy path - row_df:  {}".format(row_df))

        self.assertTrue(row_df.equals(e_row_df), (
            "\nrow_df:\n{}\ne_row_df:\n{}").format(row_df, e_row_df))

        #test that RID's are strings - using number within str
        full_df.iloc[4,0] = "13"
        full_df.iloc[5,0] = "17"
        row_df = pg.assemble_row_metadata(full_df, full_df_dims[3],
                                          full_df_dims[0], full_df_dims[2])
        logger.debug("test that RID's are strings - row_df:  {}".format(row_df))
        self.assertEqual({"13", "17"}, set(row_df.index))

    def test_assemble_col_metadata(self):
        #simple happy path
        full_df = pd.DataFrame(
            [["id", "rhd1", "rhd2", "cid1", "cid2"],
             ["chd1", "", "", "a", "b"],
             ["chd2", "", "", "50", "60"],
             ["chd3", "", "", "1.0", np.nan],
             ["rid1", "C", "D", "0.3", "0.2"],
             ["rid2", "1.0", "2.0", np.nan, "0.9"]])
        full_df_dims = [2, 2, 2, 3]
        e_col_df = pd.DataFrame([["a", 50, 1.0], ["b", 60, np.nan]],
                                index=["cid1", "cid2"],
                                columns=["chd1", "chd2", "chd3"])
        col_df = pg.assemble_col_metadata(full_df, full_df_dims[3],
                                          full_df_dims[2], full_df_dims[1])
        logger.debug("simple happy path - col_df:  {}".format(col_df))

        self.assertTrue(col_df.equals(e_col_df))

        #test that CID's are strings - using number within str
        full_df.iloc[0,3] = "13"
        full_df.iloc[0,4] = "17"
        col_df = pg.assemble_col_metadata(full_df, full_df_dims[3],
                                          full_df_dims[2], full_df_dims[1])
        logger.debug("test that CID's are strings - col_df:  {}".format(col_df))
        self.assertEqual({"13", "17"}, set(col_df.index))

    def test_assemble_data(self):
        full_df = pd.DataFrame(
            [["id", "rhd1", "rhd2", "cid1", "cid2"],
             ["chd1", "", "", "a", "b"],
             ["chd2", "", "", "55", "61"],
             ["chd3", "", "", "nah", "nope"],
             ["rid1", "C", "D", "0.3", "0.2"],
             ["rid2", "1.0", "2.0", np.nan, "0.9"]])
        full_df_dims = [2, 2, 2, 3]
        e_data_df = pd.DataFrame([[0.3, 0.2], [np.nan, 0.9]],
                                 index=["rid1", "rid2"],
                                 columns=["cid1", "cid2"], dtype = np.float32)
        data_df = pg.assemble_data(full_df, full_df_dims[3], full_df_dims[0],
                                   full_df_dims[2], full_df_dims[1])
        self.assertTrue(data_df.equals(e_data_df))

    def test_parse(self):
        # L1000 gct
        l1000_file_path = os.path.join(FUNCTIONAL_TESTS_PATH, "test_l1000.gct")
        l1000_gct = pg.parse(l1000_file_path)
        logger.debug("parse L1000 gct - l1000_gct:  {}".format(l1000_gct))
        self.assertEqual(GCToo.GCToo, type(l1000_gct))

        # Check a few values
        self.assertAlmostEqual(l1000_gct.data_df.iloc[0, 0], 11.3819, places=4,
                        msg=("The first value in the data matrix should be " +
                             "{} not {}").format("11.3819", l1000_gct.data_df.iloc[0, 0]))
        self.assertEqual(l1000_gct.col_metadata_df.iloc[0, 0], 58,
                        msg=("The first value in the column metadata should be " +
                             "{} not {}").format("58", l1000_gct.col_metadata_df.iloc[0, 0]))
        self.assertEqual(l1000_gct.row_metadata_df.iloc[0, 0], "Analyte 11",
                        msg=("The first value in the row metadata should be " +
                             "{} not {}").format("Analyte 11", l1000_gct.row_metadata_df.iloc[0, 0]))

        # P100 gct
        p100_file_path = os.path.join(FUNCTIONAL_TESTS_PATH, "test_p100.gct")
        p100_gct = pg.parse(p100_file_path)
        logger.debug("parse P100 gct - p100_gct:  {}".format(p100_gct))
        self.assertEqual(GCToo.GCToo, type(p100_gct))

        # Check a few values
        self.assertAlmostEqual(p100_gct.data_df.iloc[0, 0], 0.9182, places=4,
                        msg=("The first value in the data matrix should be " +
                             "{} not {}").format("0.9182", p100_gct.data_df.iloc[0, 0]))
        self.assertEqual(p100_gct.col_metadata_df.iloc[0, 0], "MCF7",
                        msg=("The first value in the column metadata should be " +
                             "{} not {}").format("MCF7", p100_gct.col_metadata_df.iloc[0, 0]))
        self.assertEqual(p100_gct.row_metadata_df.iloc[0, 0], 1859,
                        msg=("The first value in the row metadata should be " +
                             "{} not {}").format("1859", p100_gct.row_metadata_df.iloc[0, 0]))

        # GCT1.2
        gct_v1point2_path = os.path.join(FUNCTIONAL_TESTS_PATH, "test_v1point2_n5x10.gct")
        gct_v1point2 = pg.parse(gct_v1point2_path)
        logger.debug("parse GCT1.2 - gct_v1point2:  {}".format(gct_v1point2))
        self.assertEqual(GCToo.GCToo, type(gct_v1point2))

        # Check a few values
        self.assertAlmostEqual(
            gct_v1point2.data_df.loc["217140_s_at", "LJP005_A375_24H_X1_B19:A06"],
            6.9966, places=4)
        self.assertEqual(gct_v1point2.row_metadata_df.loc["203627_at", "Description"], "IGF1R")

        # Make sure col_metadata_df is empty
        self.assertEqual(gct_v1point2.col_metadata_df.size, 0,
                         "col_metadata_df should be empty.")

        #row_meta_only=True
        row_metadata_df = pg.parse(l1000_file_path, row_meta_only=True)
        logger.debug("row_meta_only=True - row_metadata_df.shape:  {}".format(row_metadata_df.shape))
        self.assertEqual(pd.DataFrame, type(row_metadata_df))
        self.assertEqual((978, 11), row_metadata_df.shape)

        #col_meta_only=True
        col_metadata_df = pg.parse(l1000_file_path, col_meta_only=True)
        logger.debug("col_meta_only=True - col_metadata_df.shape:  {}".format(col_metadata_df.shape))
        self.assertEqual(pd.DataFrame, type(col_metadata_df))
        self.assertEqual((377, 35), col_metadata_df.shape)

        #subsetting
        my_rids = ["218597_s_at", "214404_x_at", "209253_at"]
        my_rids_order_in_gct = ["218597_s_at", "209253_at", "214404_x_at"]
        my_cidxs = [4, 0]
        e_data_df = pd.DataFrame({"LJP005_A375_24H_X1_B19:A07": [11.04, 7.53, 6.01],
                                  "LJP005_A375_24H_X1_B19:A03": [10.45, 8.14, 4.92]},
                                 dtype=np.float32)
        e_data_df.index = pd.Index(my_rids_order_in_gct)
        e_col_meta_df = pd.DataFrame({"pert_iname": ["DMSO", "CP-724714"],
                                      "pert_id": ["DMSO", "BRD-K76908866"]}, )
        e_col_meta_df.index = pd.Index(
            ["LJP005_A375_24H_X1_B19:A03", "LJP005_A375_24H_X1_B19:A07"])

        out_g = pg.parse(l1000_file_path, rid=my_rids, cidx=my_cidxs)
        self.assertEqual(out_g.data_df.shape, (3, 2))

        # N.B. returned object should have same order as input
        pd.util.testing.assert_frame_equal(e_data_df, out_g.data_df, check_less_precise=2, check_names=False)
        pd.util.testing.assert_frame_equal(e_col_meta_df, out_g.col_metadata_df[["pert_id", "pert_iname"]], check_names=False)

    def test_parse_gct_int_ids(self):
        path = os.path.join(FUNCTIONAL_TESTS_PATH, "test_parse_gct_int_ids.gct")
        r = pg.parse(path)
        logger.debug("r:  {}".format(r))
        logger.debug("r.data_df:  {}".format(r.data_df))
        self.assertEqual({"1", "2"}, set(r.data_df.columns))
        self.assertEqual({"3", "11", "-3"}, set(r.data_df.index))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
