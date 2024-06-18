MINIWDL_CFG=./tests/miniwdl.cfg
MINIWDL_CALL_CACHE=./miniwdl_call_cache
MINIWDL_SINGULARITY_CACHE=./miniwdl_singularity_cache
MINIWDL_OUTPUT=./miniwdl_test_output

check:
	miniwdl check --strict $(wdl).wdl && echo "Workflow is valid" || echo "Workflow is invalid"

check-all:
	miniwdl check --strict *.wdl && echo "All workflows are valid" || echo "Some workflows are invalid"

check-schema:
	check-jsonschema --schemafile pbaa_parameters.schema.json tests/data/pbaa_parameters.json

list-tasks:
	miniwdl check $(wdl).wdl | grep task | sed 's/^[[:blank:]]*//' | cut -d' ' -f2

task-template:
	mkdir -p tests/inputs/templates
	miniwdl input_template --task $(task) $(wdl).wdl > tests/inputs/templates/$(wdl).$(task).template.json
	cp tests/inputs/templates/$(wdl).$(task).template.json tests/inputs/$(wdl).$(task).inputs.json

task-run:
	for i in tests/inputs/$(wdl).$(task).*inputs.json; do miniwdl run --verbose --dir $(MINIWDL_OUTPUT) --cfg $(MINIWDL_CFG) --input $$i --task $(task) $(wdl).wdl; done

template:
	mkdir -p tests/inputs/templates
	miniwdl input_template $(wdl).wdl > tests/inputs/templates/$(wdl).inputs.json
	cp tests/inputs/templates/$(wdl).inputs.json tests/inputs/$(wdl).inputs.json

run:
	miniwdl run --verbose --dir $(MINIWDL_OUTPUT) --cfg $(MINIWDL_CFG) --input tests/inputs/$(wdl).inputs.json $(wdl).wdl