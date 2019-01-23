# Must be sourced, not executed!
conda list --export > ipoolseq.reqs
conda list --explicit > ipoolseq.pkgs 
(
	echo "name: ipoolseq"
	conda env export -p "$(conda info --json | sed -n 's/^\s*"active_prefix":\s*"\(.*\)\".*$/\1/p')" \
	| egrep -v '^(name|prefix): .*$'
) > ipoolseq.yaml
