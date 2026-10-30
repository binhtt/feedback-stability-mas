.PHONY: install verify reproduce
install:
	python -m pip install -r requirements.txt
verify:
	python verify_results.py
reproduce:
	python reproduce.py
