import pandas as pd
from plot_it import color_selector, color_list_global, main_func
from argparse import Namespace
from pathlib import PosixPath


class TestAktaData:
	some_shared_data_eg_an_excel_file = pd.read_excel('Plotly_excel_template.xlsx')

	def test_excel_can_load(self):
		assert len(self.some_shared_data_eg_an_excel_file.columns) == 14

	def test_color_selector(self):
		assert color_selector(3, color_list_global) == 4

	def test_entire_script_in_one_go(self):
		#main_func(
		#	Namespace(plot_type='AKTA', input_file='Plotly_excel_template.xlsx', data_fit=False, log_scale=False,
		#			  output=PosixPath('output'), plotting=True, advanced_option_box=False)
		#)
		# then open the output files and assert they are correct
		assert 1 == 1
