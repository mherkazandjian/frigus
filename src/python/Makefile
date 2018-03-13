clean:
	rm -fvr *.pyc htmlcov

test:
	cd tests && FRIGUS_DATADIR_ROOT=../../../data PYTHONPATH=..:${PYTHONPATH} pytest -v --cov --cov-report=html .
