# Contributing to CRISPRArchitect

Thank you for your interest in contributing to CRISPRArchitect!

## How to Contribute

### Reporting Bugs
- Open a GitHub Issue with a clear description
- Include: Python version, OS, error traceback, and steps to reproduce

### Suggesting Features
- Open a GitHub Issue with the tag "enhancement"
- Describe the biological use case and expected behavior

### Submitting Code
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes with tests
4. Run the test suite: `pytest tests/ -v`
5. Submit a Pull Request

### Adding a New Nuclease
Edit `utils/constants.py` → `NUCLEASE_PARAMS` dict. Include PAM, cut type, and literature reference.

### Adding a New Cell Type
Edit `utils/constants.py` → `CELL_TYPE_PARAMS` dict. Include HDR efficiency, cell cycle, and p53 status.

### Code Style
- Type hints on all public functions
- Docstrings explaining the biology (audience: molecular biologists)
- Every parameter must cite a published source

## Running Tests

```bash
cd crisprarchitect/
pip install -e ".[dev]"
pytest tests/ -v
```

## Questions?
Open a GitHub Issue or Discussion.
