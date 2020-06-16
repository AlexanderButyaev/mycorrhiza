from .dataset import Sample
from ..exceptions import LoadingError
from .dataset import Dataset

DNA_INT = {'0': 'A'}

class Structure(Dataset):
	"""Structure formatted file reader.

		Args:
			file_path (str): Path to the file containing the data.
			diploid (bool): Dilpoid genotypes occupy 2 rows.
			n_optional_cols (int): Number of columns to ignore after the third column.
			is_str: If data is microsatellites.

	"""

	def __init__(self, file_path: str, diploid: bool = True, n_optional_cols: int = 0, is_str: bool = False):

		super().__init__(file_path, diploid, is_str)

		if n_optional_cols < 0:
			raise ValueError('Number of optional columns must be greater than zero.')
		self._n_optional_cols = n_optional_cols

		self._alpha = {}

	def _iterator(self):
		for _i, _file_path in enumerate([self._file_path] + self._extra_file_path):
			valid_func = lambda x : True
			if self._mix_indices:
				_mix_index = self._mix_indices[_i]
				if not isinstance(_mix_index, int):
					_mix_index_set = set(_mix_index)
					valid_func = lambda x : x in _mix_index_set
			with open(self._file_path) as fin:
				counter = 0
				for line in fin:

					line_a = line.strip().split()

					if self._diploid:
						line_b = fin.readline().strip().split()
						if len(line_a) != len(line_b):
							raise LoadingError('Diploid loci count mismatch for sample {0}.'.format(line_a[0]))

					geno_a = line_a[3+self._n_optional_cols:]
					geno_b = line_b[3+self._n_optional_cols:] if self._diploid else []

					if not self.is_str:
						geno_a = self._replace_missing(geno_a)
						geno_b = self._replace_missing(geno_b)

					if valid_func(counter):
						yield Sample(line_a[0], len(geno_a), population=line_a[1], known=True if line_a[2] == "1" else False), geno_a + geno_b
					counter+=1

	def _replace_missing(self, l, f='N', t='-9'):
		return [f if x == t else ['A', 'T', 'G', 'C'][int(x)] for x in l]


if __name__ == '__main__':
	pass