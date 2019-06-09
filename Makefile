.PHONY: all build install test

all:
	make build
	make install
	make test



build:
	rm -rf build/ sdist/ dist/ chemdata-*/ chemdata.egg-info/
	python setup.py sdist build
	python setup.py bdist_wheel --universal
	twine check dist/*

install:
	python setup.py install --user

travisinstall:
	python setup.py install

test:
	coverage run --source chemdata ./chemdata/test.py 
	echo `which chemdata`
	# coverage run --source chemdata `which chemdata` -h
	# coverage run --source chemdata `which chemdata` LISTSUBCOMMAND
	# coverage run --source chemdata `which chemdata` LISTSUBCOMMAND | xargs -n 1 -I [] bash -c '(coverage run --source chemdata `which chemdata` [] -h >/dev/null 2>&1 || echo ERROR: [])'
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
	cd /tmp;\
	make test'
	
upload:
	twine upload dist/*

clean:
	rm -rf venv build *.egg-info dist