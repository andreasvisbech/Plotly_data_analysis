import pandas as pd
from plot_it import color_selector, color_list_global, main_func, set_param_dict
from argparse import Namespace
from pathlib import Path


class TestAktaData:
	some_shared_data_eg_an_excel_file = pd.read_excel('Plotly_excel_template.xlsx')

	def test_excel_can_load(self):
		assert len(self.some_shared_data_eg_an_excel_file.columns) == 14

	def test_color_selector(self):
		assert color_selector(3, color_list_global) == 4

	def test_entire_script_in_one_go(self, tmp_path):
		args = Namespace(
			plot_type='akta', input_file=Path('test_data/AKTA.xlsx'), data_fit=True, log_scale=False,
			output=tmp_path, plotting=True, advanced_option_box=False)
		pdict = set_param_dict()
		main_func(args, pdict)
		# then open the output files and assert they are correct
		files = set([x.name for x in tmp_path.glob("*")])
		assert {"AKTA.tsv", "AKTA.svg", "AKTA.html"} == files
		assert len(pd.read_csv(tmp_path / "AKTA.tsv", sep='\t')) == 14530
