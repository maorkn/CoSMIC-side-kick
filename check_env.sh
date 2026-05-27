#!/bin/bash
source ~/miniforge3/bin/activate cosmic_sidekick_env
echo "=== Environment ==="
which vsearch && echo "vsearch OK" || echo "vsearch MISSING"
which barrnap && echo "barrnap OK" || echo "barrnap MISSING"
which prokka && echo "prokka OK" || echo "prokka MISSING"
which python3 && echo "python3 OK" || echo "python3 MISSING"
python3 -c "import pandas, yaml, Bio; print('ALL PYTHON DEPS OK')" 2>&1
echo "=== Done ==="
