# Required library imports
import os                      # File and directory operations
import pandas as pd           # Data manipulation and analysis
import matplotlib.pyplot as plt # Creating static visualizations
import seaborn as sns         # Statistical data visualization
import numpy as np            # Numerical operations
from typing import Tuple, List # Type hints for function signatures

class DataProcessor:
    """
    Main class for processing experimental mice data.
    Handles data loading, processing and output generation for different sample types.
    """
    
    # Dictionary mapping column keys to actual CSV column names
    COLUMNS = {
        'STRAIN': 'mouse_strain',          # Mouse strain identifier
        'EXPERIMENT': 'experiment_ID',      # Experiment identifier
        'SAMPLE_TYPE': 'sample_type',       # Type of sample collected
        'TIMEPOINT': 'timepoint',          # Experimental timepoint
        'MOUSE_ID': 'mouse_ID',            # Individual mouse identifier
        'TREATMENT': 'treatment',           # Treatment type (ABX/placebo)
        'FREQUENCY': 'frequency_live_bacteria_%',  # Percentage of live bacteria
        'EXP_DAY': 'experimental_day',      # Day relative to treatment
        'COUNTS': 'counts_live_bacteria_per_wet_g', # Absolute bacteria count
        'AGE': 'mouse_age_days',           # Mouse age in days
        'SEX': 'mouse_sex'                 # Mouse sex
    }
    
    # Available sample types in the experiment
    SAMPLE_TYPES = {
        'FECAL': 'fecal',  # Fecal samples
        'ILEAL': 'ileal',  # Ileal samples
        'CECAL': 'cecal'   # Cecal samples
    }

    def __init__(self, input_file: str):
        """
        Initialize the data processor.
        
        Args:
            input_file (str): Path to the input CSV file
        """
        self.input_file = input_file
        self.data = None  # Will hold the loaded data
        self.ensure_directories()
        
    @staticmethod
    def ensure_directories():
        """Create necessary directories for data organization."""
        for directory in ['input', 'output', 'images']:
            if not os.path.exists(directory):
                os.makedirs(directory)

    def load_data(self) -> None:
        """
        Load data from CSV file and handle potential errors.
        
        Raises:
            Exception: If file cannot be loaded
        """
        try:
            self.data = pd.read_csv(self.input_file, sep=';')
            print(f"Data successfully loaded. Dimensions: {self.data.shape}")
        except Exception as e:
            print(f"Error loading data: {e}")
            raise
        
    def process_fecal_data(self) -> pd.DataFrame:
        """
        Process fecal samples data.
        
        Returns:
            pd.DataFrame: Processed fecal data or empty DataFrame if no data found
        """
        # Filter fecal data
        fecal_data = self.data[self.data[self.COLUMNS['SAMPLE_TYPE']] == self.SAMPLE_TYPES['FECAL']].copy()
        
        if fecal_data.empty:
            print("No fecal data found")
            return pd.DataFrame()
            
        # Calculate log10 of bacteria counts
        fecal_data['log10_bacteria'] = np.log10(fecal_data[self.COLUMNS['COUNTS']])
        
        # Prepare output data
        output_data = fecal_data[[
            self.COLUMNS['MOUSE_ID'],
            self.COLUMNS['TREATMENT'],
            self.COLUMNS['EXP_DAY'],
            'log10_bacteria'
        ]].copy()
        
        # Save processed data
        output_file = 'output/fecal_data.csv'
        output_data.to_csv(output_file, index=False, sep=';')
        print(f"Fecal data saved to {output_file}")
        
        return fecal_data

    def process_organ_data(self, organ_type: str) -> pd.DataFrame:
        """
        Process organ (ileal/cecal) sample data.
        
        Args:
            organ_type (str): Type of organ data to process ('ileal' or 'cecal')
            
        Returns:
            pd.DataFrame: Processed organ data or empty DataFrame if no data found
        
        Raises:
            ValueError: If organ_type is invalid
        """
        # Validate organ type
        if organ_type not in [self.SAMPLE_TYPES['ILEAL'], self.SAMPLE_TYPES['CECAL']]:
            raise ValueError(f"Invalid organ type: {organ_type}")
            
        # Filter data for treatment end (day 0)
        organ_data = self.data[
            (self.data[self.COLUMNS['SAMPLE_TYPE']] == organ_type) &
            (self.data[self.COLUMNS['EXP_DAY']] == 0)  # Treatment end
        ].copy()
        
        if organ_data.empty:
            print(f"No data found for {organ_type} at treatment end (day 0)")
            return pd.DataFrame()
        
        # Calculate log10 of bacteria counts
        organ_data['log10_bacteria'] = np.log10(organ_data[self.COLUMNS['COUNTS']])
        
        # Prepare output data
        output_data = organ_data[[
            self.COLUMNS['TREATMENT'],
            'log10_bacteria'
        ]].copy()
        
        # Save processed data
        output_file = f'output/{organ_type}_data.csv'
        output_data.to_csv(output_file, index=False, sep=';')
        print(f"{organ_type.capitalize()} data saved to {output_file}")
        
        return organ_data

class Visualizer:
    """Class for creating data visualizations of experimental results."""
    
    def __init__(self):
        """Initialize visualization settings and style configurations."""
        plt.style.use('default')
        # Define color scheme for different treatments
        self.colors = {'ABX': 'salmon', 'placebo': 'lightblue'}
        
        # Set global plot parameters
        plt.rcParams['axes.grid'] = True
        plt.rcParams['grid.alpha'] = 0.3
        plt.rcParams['grid.linestyle'] = '--'
        plt.rcParams['font.size'] = 12

    def _calculate_axis_limits(self, data: np.ndarray, margin_percent: float = 10) -> Tuple[float, float]:
        """
        Calculate axis limits with a margin.
        
        Args:
            data: Array of values to calculate limits for
            margin_percent: Percentage of margin to add
            
        Returns:
            Tuple[float, float]: Minimum and maximum axis limits
        """
        if len(data) == 0:
            return 0, 1  # Default limits for empty data
        data_min = np.min(data)
        data_max = np.max(data)
        data_range = data_max - data_min
        margin = data_range * (margin_percent / 100)
        return data_min - margin, data_max + margin
        
    def create_fecal_plot(self, data: pd.DataFrame) -> None:
        """
        Create line plot for fecal bacteria measurements.
        
        Args:
            data: DataFrame containing fecal measurements
        """
        if data.empty:
            print("No fecal data to plot")
            return
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Configure grid appearance
        ax.grid(True, alpha=0.3, linestyle='--', color='gray')
        ax.set_axisbelow(True)
        
        # Set axis limits
        y_min, y_max = self._calculate_axis_limits(data['log10_bacteria'].values)
        ax.set_ylim(y_min, y_max)
        
        unique_days = sorted(data[DataProcessor.COLUMNS['EXP_DAY']].unique())
        ax.set_xlim(min(unique_days), max(unique_days))
        
        # Create legend entries
        for treatment, color in self.colors.items():
            plt.plot([], [], color=color, label=treatment, alpha=0.8, linewidth=1)
        
        # Plot data for each mouse
        for treatment in ['placebo', 'ABX']:
            treatment_data = data[data[DataProcessor.COLUMNS['TREATMENT']] == treatment]
            
            for mouse in treatment_data[DataProcessor.COLUMNS['MOUSE_ID']].unique():
                mouse_data = treatment_data[treatment_data[DataProcessor.COLUMNS['MOUSE_ID']] == mouse]
                plt.plot(mouse_data[DataProcessor.COLUMNS['EXP_DAY']], 
                        mouse_data['log10_bacteria'],
                        color=self.colors[treatment],
                        alpha=0.5,
                        linewidth=1)
        
        # Set labels and title
        plt.title('Fecal live bacteria', fontsize=14, pad=20)
        plt.xlabel('Washout day', fontsize=12)
        plt.ylabel('log10(live bacteria/wet g)', fontsize=12)
        
        plt.legend(frameon=True)
        plt.tight_layout()
        
        # Save plot
        plt.savefig('images/fecal_bacteria.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Fecal data plot generated")

    def create_violin_plot(self, data: pd.DataFrame, organ_type: str) -> None:
        """
        Create violin plot for organ bacteria measurements.
        
        Args:
            data: DataFrame containing organ measurements
            organ_type: Type of organ ('ileal' or 'cecal')
        """
        if data.empty:
            print(f"No {organ_type} data to plot")
            return
            
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Configure grid appearance
        ax.grid(True, alpha=0.3, linestyle='--', color='gray')
        ax.set_axisbelow(True)
        
        # Set axis limits
        y_min, y_max = self._calculate_axis_limits(data['log10_bacteria'].values)
        ax.set_ylim(y_min, y_max)
        
        # Create violin plot
        sns.violinplot(
            data=data,
            x=DataProcessor.COLUMNS['TREATMENT'],
            y='log10_bacteria',
            hue=DataProcessor.COLUMNS['TREATMENT'],
            palette=self.colors,
            inner='box',
            legend=False
        )
        
        # Add individual data points
        sns.stripplot(
            data=data,
            x=DataProcessor.COLUMNS['TREATMENT'],
            y='log10_bacteria',
            color='gray',
            alpha=0.6,
            size=5,
            jitter=0.2
        )
        
        # Set labels and title
        plt.title(f'{organ_type.capitalize()} live bacteria', fontsize=14, pad=20)
        plt.xlabel('Treatment', fontsize=12)
        plt.ylabel('log10(live bacteria/wet g)', fontsize=12)
        
        # Add legend
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in self.colors.values()]
        labels = list(self.colors.keys())
        plt.legend(handles, labels, loc='lower right')
        
        plt.tight_layout()
        plt.savefig(f'images/{organ_type}_bacteria.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"{organ_type.capitalize()} data plot generated")

def main():
    """
    Main execution function.
    Handles file processing and visualization pipeline.
    """
    try:
        # Check for input files
        csv_files = [f for f in os.listdir('input') if f.endswith('.csv')]
        
        # Validate input file count
        if len(csv_files) == 0:
            raise FileNotFoundError("No CSV file found in input directory")
        elif len(csv_files) > 1:
            raise ValueError("Please put only one CSV file in input directory")
            
        # Process input file
        input_file = os.path.join('input', csv_files[0])
        print(f"Processing file: {csv_files[0]}")
        
        # Initialize processor and load data
        processor = DataProcessor(input_file)
        processor.load_data()
        
        # Create visualizations
        visualizer = Visualizer()
        
        # Process and plot fecal data
        fecal_data = processor.process_fecal_data()
        if not fecal_data.empty:
            visualizer.create_fecal_plot(fecal_data)
        
        # Process and plot organ data
        for organ in [DataProcessor.SAMPLE_TYPES['CECAL'], DataProcessor.SAMPLE_TYPES['ILEAL']]:
            organ_data = processor.process_organ_data(organ)
            if not organ_data.empty:
                visualizer.create_violin_plot(organ_data, organ)
            
        print("Processing completed successfully!")
            
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()