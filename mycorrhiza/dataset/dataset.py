from ..exceptions import LoadingError
import numpy as np

class Sample:

	def __init__(self, identifier: str, num_loci: int, population=None, known: bool = True):
		self._identifier = str(identifier)
		self._population = population
		self._known = known
		self._num_loci = num_loci

	@property
	def identifier(self):
		return self._identifier

	@property
	def num_loci(self):
		return self._num_loci

	@property
	def population(self):
		return self._population

	@property
	def flag(self):
		return self._known


class Dataset:

	def __init__(self, file_path, diploid: bool = True, is_str: bool = True):

		self._file_path = file_path
		self._diploid = diploid
		self._is_str = is_str
		self._samples = []
		self._num_loci = None
		self._extra_file_path = []
		self._mix_indices = []

	def _iterator(self):
		return []

	@property
	def iterator(self):
		return self._iterator

	def _add_sample(self, sample: Sample):

		if self._num_loci is None:
			self._num_loci = sample.num_loci
		elif self._num_loci != sample.num_loci and self._num_loci is not None:
			raise LoadingError('Mismatch in the number of loci.')

		self._samples.append(sample)

	def load(self):
		self._samples = []

		for sample, genotype in self._iterator():
			self._add_sample(sample)

		self._mix_indices.append( len(self._samples) )
		self.statistics()

	def statistics(self):
		print('Loaded {0} samples with {1} loci.'.format(self.num_samples, self.num_loci))


	@property
	def num_samples(self):
		return len(self._samples)

	@property
	def num_loci(self):
		return self._num_loci

	@property
	def haploid(self):
		return not self._diploid

	@property
	def diploid(self):
		return self._diploid

	@property
	def is_str(self):
		return self._is_str

	@property
	def populations(self):
		return [sample.population for sample in self._samples]

	@property
	def identifiers(self):
		return [sample.identifier for sample in self._samples]

	@property
	def flags(self):
		return [sample.flag for sample in self._samples]

	def setSamplePopulation(self, index, population):
		self._samples[index]._population = population

	def _microsatellite_distances(self):
		matrix = np.array([x[1] for x in self._iterator()])
		distances = np.zeros([matrix.shape[0]]*2)
		for i in range(matrix.shape[0]):
			for j in range(i+1, matrix.shape[0]):
				mtx = matrix[[i,j],:]
				mtx = mtx[:, ~np.isin(mtx[0,:],['-9','000']) & ~np.isin(mtx[1,:],['-9','000'])]
				sum_vals = (mtx[0,:] != mtx[1,:]).sum()
				len_vals = mtx.shape[1]
				distances[i, j] = sum_vals/len_vals if len_vals != 0 else 0
		distances = np.maximum(distance, distances.T)
		return distances

	def concatenate(self, *others):
		# check compatibility before loading
		# @TODO: self._num_loci - align with original file...
		for _dataset in others:
			if not (self._diploid == _dataset.diploid and
					self._is_str == _dataset._is_str and
					self._num_loci == _dataset._num_loci):
				raise ValueError(
					"Original dataset and appendix must have the same properties." +
					" (Source: %s)" % _dataset._file_path
				)
		# mix with already existing identifiers
		other_indices = []
		sample_no = original_sample_size = self.num_samples
		for _dataset in others:
			indices = []
			identifiers = self.identifiers
			are_distinct = True
			sub_mix_indices = []
			for _i, _sample in enumerate(_dataset._samples):
				try:
					indx = identifiers.index(_sample.identifier)
					are_distinct = False
				except ValueError:
					indx = int(sample_no)
					self._samples.append(_sample)
					sub_mix_indices.append(_i)
					sample_no += 1

				indices.append(indx)
			if are_distinct:
				self._mix_indices.append(_dataset.num_samples)
			else:
				self._mix_indices.append(sub_mix_indices)
			self._extra_file_path.append(_dataset._file_path)
			other_indices.append(indices)
		self.statistics()
		return original_sample_size, other_indices
