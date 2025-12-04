# Contributing to Lysozyme Annotation Pipeline

Thank you for considering contributing to this project!

## How to Contribute

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/amazing-feature`)
3. **Make your changes**
4. **Test thoroughly** with example data
5. **Commit your changes** (`git commit -m 'Add amazing feature'`)
6. **Push to the branch** (`git push origin feature/amazing-feature`)
7. **Open a Pull Request**

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/lysozyme_pipeline.git
cd lysozyme_pipeline

# Install dependencies
pip install -r requirements.txt

# Run tests
python pipeline.py -g examples/example_genome.fasta -l examples/reference_lysozymes.fasta -o test_output/
```

## Code Style

- Follow PEP 8 conventions
- Add docstrings to all functions
- Use type hints where appropriate
- Keep functions focused and modular

## Testing

Before submitting:
- Test with single genome mode
- Test with batch mode
- Verify output files are generated correctly
- Check for errors in logs

## Reporting Bugs

Please include:
- Pipeline version
- Python version
- Operating system
- Complete error message
- Minimal example to reproduce

## Suggesting Enhancements

Open an issue describing:
- Use case
- Expected behavior
- Potential implementation approach
