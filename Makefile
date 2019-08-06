.PHONY: all reqs build install test

pes_parent_dir:=$(shell pwd)/$(lastword $(MAKEFILE_LIST))
pes_parent_dir:=$(shell dirname $(pes_parent_dir))

Project=$(shell basename $(pes_parent_dir))

all:
	make reqs
	make build
	make install
	make test

reqs:
	pipreqs --help >/dev/null 2>&1 || pip3 install pipreqs || pip3 install pipreqs --user
	pipreqs --force $(Project)
	mv $(Project)/requirements.txt .
	sed -i 's/==/>=/g' requirements.txt
	sed -i 's/numpy.*/numpy/g' requirements.txt
	sed -i 's/psutil.*/psutil/g' requirements.txt
	cat requirements.txt 

build:
	rm -rf build/ sdist/ dist/ $(Project)-*/ $(Project).egg-info/
	# python setup.py sdist build
	python setup.py bdist_wheel --universal
	twine check dist/*

install:
	cd /tmp; pip uninstall -yy $(Project); cd -; python setup.py install || python setup.py install --user

test:
	bash -c "export PYTHONPATH="$(PYTHONPATH):$(PWD)"; coverage run --source $(Project) ./tests/test.py" 
	echo `which $(Project)`
	# coverage run --source $(Project) `which $(Project)` -h
	# coverage run --source $(Project) `which $(Project)` LISTSUBCOMMAND
	# coverage run --source $(Project) `which $(Project)` LISTSUBCOMMAND | xargs -n 1 -I [] bash -c '(coverage run --source $(Project) `which $(Project)` [] -h >/dev/null 2>&1 || echo ERROR: [])'
	coverage report -m

test_env:
	bash -c ' \
	rm -rf venv; \
	virtualenv venv; \
	source venv/bin/activate; \
	which python; \
	python --version; \
	pip install -r requirements.txt; \
	make build; \
	make travisinstall; \
	make test'
	
upload:
	twine upload dist/*

clean:
	rm -rf venv build *.egg-info dist
