clean:
	rm -fvr *.pyc htmlcov

test:
	cd tests && FRIGUS_DATADIR_ROOT=../../../data PYTHONPATH=..:${PYTHONPATH} py.test -v --cov --cov-report=html .
