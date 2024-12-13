# Projet-CY-Antibio-Tech

Python program for analyzing and visualizing the impact of antibiotics on mice intestinal microbiota.

## Setup and Installation

### Prerequisites
- Python 3.x
- Required packages: pandas, numpy, matplotlib, seaborn, scipy

### Installation
1. For macOS:
```bash
pip3 install pandas numpy matplotlib seaborn scipy
```

2. For Windows:
```bash
pip install pandas numpy matplotlib seaborn scipy
```

## Usage Instructions

1. Put your CSV data file in the `input` folder (only one file at a time)
2. Run the program:
   - macOS: `python3 antibio_analysis.py`
   - Windows: `python antibio_analysis.py`
3. Check results in:
   - `output/` for processed CSV files
   - `images/` for generated graphs

## Data Processing Details

The program processes three types of data:
1. Fecal data: tracked throughout the experiment
2. Cecal data: measured at treatment end (day 0)
3. Ileal data: measured at treatment end (day 0)

## Project Limitations

### Known Limitations
- Only processes one CSV file at a time
- CSV must use semicolon (;) as separator
- No handling of missing values
- Fixed color scheme for graphs (salmon for ABX, lightblue for placebo)
- Graph labels are in English only
- Organ data (cecal and ileal) are only processed for day 0 (treatment_end)

### Successfully Implemented
- Automatic file detection in input directory
- Processing of fecal, ileal, and cecal data
- Generation of line plots for fecal data
- Generation of violin plots for organ data
- Support for variable number of mice
- Automatic axis scaling
- Empty data handling

### Potential Improvements
- Add data validation
- Add support for multiple languages
- Implement additional plot types
- Add statistical analysis
- Add configuration options
- Add support for different experimental days for organ data

Note: This program assumes your CSV file follows the specified structure with correct column names and data types.
