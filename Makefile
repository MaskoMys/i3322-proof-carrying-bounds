export PYTHONDONTWRITEBYTECODE := 1

.PHONY: verify verify-all manifest refresh-manifest primary second mutations optimized optimized-pristine optimized-corruption paper clean

verify: manifest primary

verify-all: verify second mutations optimized

manifest:
	python3 tools/verify_manifest.py

refresh-manifest:
	python3 tools/build_manifest.py

primary:
	python3 verifier/verify_all.py

second:
	python3 verifier/verify_independent.py

mutations:
	python3 tests/test_mutations.py

optimized: optimized-pristine optimized-corruption

optimized-pristine:
	bash tests/test_optimized_pristine.sh

optimized-corruption:
	bash tests/test_optimized_corruption.sh

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode i3322_exact_reproducible_bounds.tex

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type f \( -name '*.pyc' -o -name '*.aux' -o -name '*.log' -o -name '*.out' -o -name '*.fls' -o -name '*.fdb_latexmk' \) -delete
