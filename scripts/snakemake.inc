import os
import fnmatch
import traceback
from snakemake.io import contains_wildcard, regex

if 'TMPDIR' in os.environ and 'USER' in os.environ:
	SCRATCH_DIR = os.path.join(os.environ['TMPDIR'], os.environ['USER'])
elif 'TMPDIR' in os.environ:
	SCRATCH_DIR = os.environ['TMPDIR']
else:
	SCRATCH_DIR = '/tmp'

def strip_ext(path, ext):
	if not path.endswith(ext):
		return path
	else:
		return path[:len(path)-len(ext)]

def get_options(rule, input=None, files=None, wildcards=None, required=False):
	try:
		if input is not None:
			if files is not None or wildcards is not None:
				raise ValueError('get_options: Either input, or files and wildcards must be specified')
			files = [file for file in input]
		elif files is not None and wildcards is not None:
			files = [ files ] if isinstance(files, str) else files
			if input is not None:
				raise ValueError('get_options: Either input, or files and wildcards must be specified')
			files = [file.format(**wildcards) for file in files]
		if rule not in config.keys():
			raise RuntimeError("CFGFILE: No configuration entry for rule %s" % rule)
		if not isinstance(config[rule], list):
			raise RuntimeError("CFGFILE: Options for rule %s consist of a %s instead of a list of dictionaries" %
			                   (rule, type(e)))
		for e in config[rule]:
			if not isinstance(e, dict):
				raise RuntimeError("CFGFILE: Option for rule %s consist of a %s instead of a dictionary" %
						   (rule, type(e)))
			if 'file' not in e or 'opts' not in e:
				raise RuntimeError("CFGFILE: Option for rule %s is missing a 'file' or 'opts' specification" %
						   (rule))
			for file in files:
				if 'file' not in e.keys():
					raise RuntimeError("CFGFILE: Missing 'files' specification for rule %s" % rule)
				if contains_wildcard(e['file']):
					wildcards = glob_wildcards(e['file'], files=[file])
					if len(wildcards[0]) == 0:
						continue
				else:
					wildcards = None
					if not fnmatch.fnmatch(file, e['file']):
						continue
				if 'opts' not in e.keys():
					raise RuntimeError("CFGFILE: Missing 'opts' for files '%' and rule %s" %
					                   (e['file'], rule))
				if wildcards is None:
					return e['opts']
				else:
					return expand(e['opts'], **wildcards._asdict())[0]
		if required:
			raise RuntimeError("CFGFILE: No matching 'files' specifiction for rule %s and file %s" % (rule, file))
		else:
			return ''
	except Exception as e:
		raise
#		print("WARNING: get_options for %s failed with: %s" % (rule, str(e)), file=sys.stderr)
#		traceback.print_exc(file=sys.stderr)
#		sys.stderr.flush()


def config_options(rule, required=False):
	return lambda wildcards, input: get_options(rule, input=input, required=required)

def config_input_file(rule, related_inputs, stripext="", ext=""):
	return lambda wildcards : strip_ext(get_options(rule, files=related_inputs, wildcards=wildcards, required=True), stripext) + ext

def get_scratch(output):
	if len(output) == 1:
		ofile = str(output)
	else:
		ofile = next(output.items().__iter__())[1]
	return os.path.join(SCRATCH_DIR, '%s@%d' % (ofile.replace('/', '-'), os.getpid()))

default_scratch = lambda wildcards, output: get_scratch(output)

