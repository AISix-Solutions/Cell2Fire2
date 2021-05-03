"""
Before you run this, run

  python main.py --input-instance-folder ../data/Sub40x40/ --output-folder ../results/Sub40x40 --sim-years 2 --nsims 10 --finalGrid --weather random --nweathers 100 --Fire-Period-Length 1.0 --ROS-CV 0.0 --seed 123 --IgnitionRad 0 --stats --output-messages --ROS-Threshold 0 --HFI-Threshold 0 --heuristic 3 --customValue="../data/Sub40x40/values40x40.csv"

That will create the input files. Then we run basically the same command under the tester, but
with --onlyProcessing to avoid calling the compiled C++ code.
# TBD: this test module should take the command line as an arg, so it does not have to
#      be pasted into the github action yml file and here. It should just be in the action
#      yml file and read here, with --onlyProcessing appended here.

"""
import unittest
import os.path
import datetime
from cell2fire.utils.ParseInputs import make_parser
from cell2fire.Cell2FireC_class import *
import cell2fire  # for path finding
import pandas as pd


p = str(cell2fire.__path__)
l = p.find("'")
r = p.find("'", l+1)
cell2fire_path = p[l+1:r]
data_path = os.path.join(cell2fire_path, "..","data")

def _readme_list():
    datadir = os.path.abspath(os.path.join(data_path, "Sub40x40"))
    resultsdir = os.path.abspath(os.path.join("..", "results", "Sub40x40"))
    customValueDir = os.path.abspath(os.path.join(datadir,  "values40x40.csv"))
    baselist = ["--input-instance-folder", datadir,
                "--output-folder", resultsdir,
                "--sim-years", "2",
                "--nsims",  "10",
                "--finalGrid",
                "--weather", "random",
                "--nweathers", "100",
                "--Fire-Period-Length", "1.0",
                "--ROS-CV", "0.0",
                "--seed", "123",
                "--IgnitionRad", "0",
                "--stats",
                "--output-messages",
                "--ROS-Threshold", "0",
                "--HFI-Threshold", "0",
                "--heuristic", "3",
                "--customValue", "../data/Sub40x40/values40x40.csv"]
    return baselist


class TestMain(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        pass

    def test_readme_cmd(self):
        print("Running ", str(self.id()).split('.')[2])
        parser = make_parser()
        cmdlist = _readme_list()
        cmdlist.append("--onlyProcessing")
        # only-processing means something else has to make the data!!! (no call to C++)

        args = parser.parse_args(cmdlist)
        env = Cell2FireC(args)  # see main.py
        env.stats()

        results_path = os.path.abspath(os.path.join("..", "results", "Sub40x40"))
        csv_path = os.path.join(results_path, "Stats", "FinalStats.csv")
        df = pd.read_csv(csv_path)

        self.assertAlmostEqual(df['Burned'][0], 103), "TEST ERROR"
        self.assertAlmostEqual(df['Burned'][6], 570), "TEST ERROR"
        self.assertAlmostEqual(df['Burned'][9], 80), "TEST ERROR"


if __name__ == "__main__":

    """
    # first: full run to create data files for testing  
    Unfortunately, we can run't under os.system because the call to compiled code fails....
    cmdstr = "python main.py "
    for es in _readme_list():
        cmdstr += es+' '
    print("cmdstr=",cmdstr)
    os.system(cmdstr)
    """
    # something else needs to create the files
    # second: do the tests
    unittest.main()
