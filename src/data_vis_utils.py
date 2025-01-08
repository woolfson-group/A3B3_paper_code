import json
import math
import re
from typing import Union, Tuple, Dict, List
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from PIL import Image
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error

def read_hplc_csv(sample_file: Path, blank_file: Path, peptide: str, wavelength: int) -> pd.DataFrame:
    """
    Reads HPLC data from a sample and a blank CSV file, processes the data to 
    calculate baseline-corrected and normalized intensity values, and returns 
    the processed data as a pandas DataFrame.
    """

    # Read the sample CSV file into a DataFrame
    df_sample = pd.read_csv(sample_file, names=['t_min', 'int'], skiprows=1)

    # Read the blank CSV file into a DataFrame
    df_blank = pd.read_csv(blank_file, names=['t_min', 'int'], skiprows=1)

    # Add peptide and wavelength information as constant columns
    df_sample['peptide'] = peptide
    df_sample['wavelength'] = wavelength

    # Calculate baseline-corrected intensity (sample minus blank)
    df_sample['int_zero'] = df_sample['int'] - df_blank['int'].values

    # Calculate intensity with the minimum sample intensity subtracted (no blank correction)
    df_sample['int_zero_no_blank'] = df_sample['int'] - df_sample['int'].min()

    # Calculate min-max normalized baseline-corrected intensity
    df_sample['int_min_max'] = (
        (df_sample['int_zero'] - df_sample['int_zero'][0]) / 
        (df_sample['int_zero'].max() - df_sample['int_zero'][0])
    )

    return df_sample


class CircularDiochroismFunctions:
    """
    A series of functions to tidy, transform and visualise CD data from JASCO 810/815
    spectrophotometers.
    """

    @staticmethod
    def calc_mre(deg: float, conc: float, amides: int, pathlength_cm: float) -> float:
        """
        Calculate the mean residue ellipticity (MRE) from ellipticity in degrees.

        Formulation taken from:
        https://www.photophysics.com/faqs/methods-techniques/cd-units-faqs/

        Args:
            deg (float): The ellipticity (degrees).
            conc (float): Concentration of peptides (molar).
            amides (int): Number of amide bonds (residues - 1).
            pathlength_cm (float): The length of the transmission path (cm).

        Returns:
            float: The calculated MRE in units of deg cm² dmol⁻¹ res⁻¹, or 
            np.nan if inputs are invalid.
        """
        # Return np.nan if any of the critical parameters are non-positive
        if conc <= 0 or amides <= 0 or pathlength_cm <= 0:
            return np.nan

        # Calculate MRE
        return 100 * (deg / (conc * amides * pathlength_cm))

    @staticmethod
    def calc_mre_row(row):
        """Define a wrapper function for row-wise application converting units"""
        return CircularDiochroismFunctions.calc_mre(
            deg=row['abs_zero_mdeg'] / 1_000,
            conc=row['conc_uM'] / 1_000_000,
            amides=row['amides'],
            pathlength_cm=row['path_mm'] / 10,
        )

    @staticmethod
    def calc_frac_helix(mre222: float) -> float:
        """
        Calculate the fractional helicity based on MRE at 222 nm.

        Equation taken from:
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4395514/

        Args:
            mre222 (float): Mean residue ellipticity at 222 nm (deg cm² dmol⁻¹).

        Returns:
            float: The calculated fractional helicity (0 to 1).
        """

        mre_ref_min = -39500  # MRE for fully helical
        mre_ref_max = -3000   # MRE for fully random coil

        return (mre222 - mre_ref_max) / (mre_ref_min - mre_ref_max)

    @staticmethod
    def extract_cd_info(file_path: Union[str, Path]) -> Dict[str, Union[str, float]]:
        """
        Extracts and returns relevant information from a file name in the given file path.
        
        The file names should follow this format:
            {date}_{type}_{peptide_info}_{amides}amides_{path_mm}mm_{temp_C}C_{solvent}.txt

        Args:
            file_path (Union[str, Path]): The file path containing the file to be processed.

        Returns:
            dict: A dictionary containing the extracted information.

        Raises:
            ValueError: If the file name format is incorrect.
        """

        # Convert file_path to Path if it's a string
        file_path = Path(file_path) if isinstance(file_path, str) else file_path

        file_name = file_path.name

        # Define regex patterns
        patterns = {
            "date": r"^(\d+)",
            "type": r"^\d+_([a-z]+)",
            "amides": r"(\d+)amides",
            "path_mm": r"(\d+)mm",
            "temp_range_C": r"(\d+-*\d*)C",
            "solvent": r"\d+C_(.*).txt",
            "peptide_info": r"^\d+_[a-z]+_(.*)_\d+amides",
        }

        # Extract basic information
        try:
            info_dict = {
                "date": re.search(patterns["date"], file_name).group(1),
                "type": re.search(patterns["type"], file_name).group(1),
                "amides": float(re.search(patterns["amides"], file_name).group(1)),
                "path_mm": float(re.search(patterns["path_mm"], file_name).group(1)),
                "temp_range_C": re.search(patterns["temp_range_C"], file_name).group(1),
                "solvent": re.search(patterns["solvent"], file_name).group(1).replace("_", " "),
            }
        except AttributeError:
            raise ValueError(f"File name format is incorrect: {file_name}")

        # Handle "blank" type files
        if info_dict["type"] == "blank":
            info_dict["related_id"] = "water"
            info_dict["conc_uM"] = 0
        else:
            # Extract peptide information
            peptide_info = re.search(patterns["peptide_info"], file_name).group(1)
            info_split = peptide_info.split("_")

            peptides = [info_split[i] for i in range(1, len(info_split), 2)]
            amide_concs = [float(info_split[i][:-2]) for i in range(0, len(info_split), 2)]

            info_dict["related_id"] = "+".join(peptides)
            info_dict["conc_uM"] = sum(amide_concs)

            # Add individual concentrations
            for i, conc in enumerate(amide_concs):
                key = f"conc_uM_{i+1}"
                info_dict[key] = conc if not math.isnan(conc) else 0

        return info_dict

    @staticmethod
    def subtract_blank(row: pd.Series, blank_data: pd.DataFrame) -> float:
        """
        This function adjusts the 'abs_mdeg' value of a given row by subtracting the 
        'abs_mdeg' value of the matching blank row.

        Args:

            row : pd.Series
                The row of data for which the 'abs_mdeg' value needs to be adjusted. 
                This row must have 'type', 'date', 'solvent', 'temp_range_C', and 
                'path_mm' fields.

            blank_data : pd.DataFrame
                The DataFrame containing blank data, which must have 'date', 'solvent',
                'temp_range_C', 'path_mm', and 'abs_mdeg' fields.

        Returns:
            float: The adjusted 'abs_mdeg' value for the given row.
        """

        # Define the list of types that require blank subtraction
        types_requiring_blank_subtraction = ['pre', 'post']

        # Check if the type of the row is in the list
        if row['type'] in types_requiring_blank_subtraction:
            # Find the matching blank by comparing 'date', 'solvent', 'temp_range_C', and 'path_mm'
            matching_blank = blank_data[
                (blank_data['date'] == row['date']) &
                (blank_data['solvent'] == row['solvent']) &
                (blank_data['temp_range_C'] == row['temp_range_C']) &
                (blank_data['path_mm'] == row['path_mm'])
            ]
            # If there's a matching blank, adjust the 'abs_mdeg' value
            if not matching_blank.empty:
                blank_value = matching_blank.iloc[0]['abs_mdeg']
                return row['abs_mdeg'] - blank_value

        # If the row's type is not in the list or no matching blank was found, return the original 'abs_mdeg' value
        return row['abs_mdeg']

    @staticmethod
    def make_cd_df(file_path: Union[str, Path]) -> pd.DataFrame:

        info = CircularDiochroismFunctions.extract_cd_info(file_path)

        df = pd.read_csv(
            file_path,
            sep='\t',
            skiprows=19,
            header=None,
            names=['wavelength_nm', 'abs_mdeg', 'ht_V']
        )

        # Change column header to "T_C" if it's a melt or cool
        if info['type'] in ["melt","cool"]: 
            df.columns = ['T_C', 'abs_mdeg', 'ht_V']

        # Append values from the Series to the DataFrame
        df = df.assign(**info)

        return df
    
    @staticmethod
    def determine_pep_type(s:str):
        if '+' in s:
            return 'AB'
        if 'Hex2-B' in s:
            return 'B' 
        if 'Hex2-A' in s:
            return 'A'
        return 'water'
    
    @staticmethod
    def determine_cnf_type(s:str):
        if 'n4CF' in s:
            return 'n4CF'
        if 'c4CF' in s:
            return 'c4CF'
        return 'NA'

class AUCanalysis:

    @staticmethod
    def parse_sv(data_path: Path, mono_mass: int) -> Tuple[pd.DataFrame, Image.Image, str, dict]:

        # Create a dictionary of all files in the directory, keyed by their suffixes
        files = {f.suffix: f for f in data_path.iterdir() if f.is_file()}

        # Open the info.txt file and read its content
        with files['.txt'].open() as info_file:
            info = info_file.read().replace('\n', ' ')

        # Extract info values from info string using regex patterns
        info_values = [float(re.search(pattern, info).group(1)) for pattern in ['vbar = (\d.\d*)', 'Mw = (\d*) Da ', 'sw = (\d.\d{3}) S', 'sw\W20,w\W = (\d.\d{3}) S', 'ratio = (\d.\d{3})']]

        # Assign values to variables
        vbar, mass, sw, sw20w, ff0 = info_values

        # Calculate oligo value
        oligo = round(mass/mono_mass, 1)

        # Load data from csv file and assign related_id
        data_df = pd.read_csv(files['.csv'], names=['Sedimentation Coefficient (S)', 'c(s)'])

        # Prepare parameters for return
        params = {
            'vbar':vbar, 
            'mass':mass, 
            'oligo':oligo, 
            'sw':sw, 
            'sw20w':sw20w, 
            'ff0':ff0, 
            'mono_mass':mono_mass
        }

        # Open image file
        bitmap = Image.open(files['.png'])

        # Return dataframe, bitmap, info string and parameters dictionary
        return data_df, bitmap, params

    def header_csv(csv_path, num_channels, alt=False):
        
        # get the csv and get values
        df_csv = pd.read_csv(csv_path, header=None)
        if df_csv.iloc[0,0] != 'channel_1':
            # generate file path and directory
            directory = csv_path.parent

            # get the different speeds
            RA_files = sorted(directory.rglob('0*.RA*'))
            
            speeds = []
            for p in RA_files[::2]:
                with p.open() as f:
                    speed = int(re.search('R \d+ \d+\.\d (\d*)', f.read()).group(1))
                    speeds.append(speed)

            # generate the correct multiindex for the new column headers
            tuples = []
            if alt:
                for data_type in ['raw', 'fit', 'res']:
                    for channel in range(1,num_channels+1):
                        for speed in speeds:
                            for xy in ['radius', 'abs']:
                                tuples.append((f'channel_{channel}', data_type, speed, xy))
            else:
                for channel in range(1,num_channels+1):
                    for data_type in ['raw', 'fit', 'res']:
                        for speed in speeds:
                            for xy in ['radius', 'abs']:
                                tuples.append((f'channel_{channel}', data_type, speed, xy))

            # add columns to csv and save
            df_csv.columns = pd.MultiIndex.from_tuples(tuples)
            df_csv.to_csv(csv_path, index=False)
        else:
            print('File already has formatted headers.')
            
    def se(p: Path, mono_mass: float):
        """
        outpath: output file path without the suffix.    
        """
        # function parameters
        dir_path=p.parent

        # get mc iterations mass mean 
        mc_df = pd.read_csv(dir_path/'mc.dat-raw', header=None, delimiter='\t')
        mass = int(mc_df.iloc[:,0].mean())

        # get and tidy csv file
        df = pd.read_csv(p, header=[0,1,2,3])
        df = df.replace(0, np.nan)
        for col in df.iloc[:,range(0,df.shape[1],2)].columns:
            df.loc[:,col] = df.loc[:,col]**2-min(df.loc[:,col]**2) # converts radius to R-R^2

        df_stack = df.stack([0,1,2])
        df_tidy = df_stack.reset_index().iloc[:,1:].dropna()
        df_tidy.columns = ['channel', 'data_type', 'speed', 'radius', 'ab']
        
        df_tidy.to_csv(p.parent / f'{p.stem}_tidy.csv', index=False)
        
        results_dict = {
            'data_set': p.name,
            'mono_mass_Da':mono_mass,
            'molecular_mas_Da':mass,
            'oligomeric_state': round(mass/mono_mass, 2)
        }
        
        with open(p.parent / f'{p.stem}_tidy.json', 'w') as f:
            json.dump(results_dict, f, indent=4)

class DPHBindingAnalysis:

    def __init__(self, file_path: Path):
        """
        Initialise the DPHBindingAnalysis object with an Excel file path.
        
        Parameters:
        - file_path (str): The path to the Excel file containing the data.
        """
        self.file_path = file_path
        self.data = pd.DataFrame()
        self.filtered_data = pd.DataFrame()
        self.normalised_data = pd.DataFrame()
        self.grouped_data = pd.DataFrame()
        self.fit_metrics = {}

    def load_binding_data(self) -> pd.DataFrame:
        """Load the Excel file into a DataFrame."""
        return pd.read_excel(
            self.file_path,
            skiprows=11,
            usecols=[3, 4],
            names=['int_AU', 'conc_uM']
        )

    def filter_by_iqr(self, group: pd.DataFrame, column: str) -> pd.DataFrame:
        """
        Filter rows outside the IQR for a specific column in a DataFrame group.
        
        Parameters:
        - group (pd.DataFrame): DataFrame containing the group data.
        - column (str): The name of the column to apply IQR filtering on.
        
        Returns:
        - pd.DataFrame: The filtered DataFrame.
        """
        Q1 = group[column].quantile(0.25)
        Q3 = group[column].quantile(0.75)
        IQR = Q3 - Q1
        return group.query(f'(@Q1 - 1.5 * @IQR) <= {column} <= (@Q3 + 1.5 * @IQR)')

    def min_max_normalise(self, df: pd.DataFrame, column: str) -> pd.DataFrame:
        """Apply min-max normalisation to a specific column in the DataFrame."""
        min_val = df[column].min()
        max_val = df[column].max()
        df[f'norm_{column}'] = (df[column] - min_val) / (max_val - min_val)
        return df

    def tidy_data(self):
        """Main function to load, filter, and normalise data."""
        self.data = self.load_binding_data()
        self.filtered_data = self.data.groupby('conc_uM').apply(self.filter_by_iqr, column='int_AU').reset_index(drop=True)
        self.normalised_data = self.min_max_normalise(self.filtered_data, 'int_AU')
    
    @staticmethod
    def single_site(x: np.ndarray, bmax: float, Kd: float) -> np.ndarray:
        """Model function for a single site binding."""
        return bmax * x / (Kd + x)

    def fit_and_evaluate_model(self):
        """
        Fit the model to the data and compute evaluation metrics.
        Update the fit_metrics attribute.
        """
        self.grouped_data = self.normalised_data.groupby('conc_uM')['norm_int_AU'].agg(['mean', 'sem']).reset_index()
 
        params, params_err = curve_fit(self.single_site, self.grouped_data['conc_uM'], self.grouped_data['mean'])
        fitted_values = self.single_site(self.grouped_data['conc_uM'], *params)
        
        rmse = np.sqrt(mean_squared_error(self.grouped_data['mean'], fitted_values))
        residuals = self.grouped_data['mean'] - fitted_values
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.grouped_data['mean'] - np.mean(self.grouped_data['mean'])) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        
        self.fit_metrics = {'params': params, 'params_err':params_err, 'RMSE': rmse, 'R^2': r_squared}
        
    def plot_data_and_fit(self, outfile:str, params: dict):
        """Plot the data and the model fit."""
        
        sns.set(**params)
        x_values = np.linspace(0, max(self.grouped_data['conc_uM']), 151)
        y_values = self.single_site(x_values, *self.fit_metrics['params'])
        plt.errorbar(self.grouped_data['conc_uM'], self.grouped_data['mean'], yerr=self.grouped_data['sem'], fmt='o', label='Data', color='k', mec='white', ms=7)
        plt.plot(x_values, y_values, label='Fit', linewidth=2, color='k')
        plt.xlabel('Concentration (uM)')
        plt.ylabel('Normalised Intensity (AU)')
        
        plt.tight_layout()
        plt.savefig(outfile, dpi=400, transparent=True)
        plt.show()
        
    def write_log(self) -> None:
        """Write analysis summary to a JSON log file."""
        
        # Create log data
        log_data = {
            'date_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'file_path': str(self.file_path),
            'num_data_points': len(self.data),
            'num_data_points_after_filter': len(self.normalised_data),
            'fit_metrics': {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in self.fit_metrics.items()},
        }
        
        # Write processed data to CSV file
        self.normalised_data.to_csv(self.file_path.with_suffix('.csv'), index=False)
        
        # Write log data to JSON file
        with open(self.file_path.with_suffix('.json'), 'w') as f:
            json.dump(log_data, f, indent=4)

class AF3Analysis:
    """A class of functions for processing AF3 data tidying / visualisation."""

    @staticmethod
    def parse_AF3_scores(json_path:Path, cif_path:Path) -> Tuple[List[List], List]:

        # load score json file
        with open(json_path, 'r') as json_obj:
            json_data = json.load(json_obj)

        # get CA plddts from cif files
        ca_plddts = []
        with open(cif_path, 'r') as cif_obj:
            for line in cif_obj.readlines():
                if line.startswith('ATOM') and ' CA ' in line:
                    ca_plddt = float(line.split()[-4])
                    ca_plddts.append(ca_plddt)

        return json_data['pae'], ca_plddts

    @staticmethod
    def AF3_plddt_plot(plddt_data:List):
        plddt_ax = sns.lineplot(plddt_data, color='k', lw=2, zorder=1)

        # Add coloured patches
        plddt_ax.axhspan(0, 50, facecolor='#FFA500', zorder=0)   
        plddt_ax.axhspan(50, 70, facecolor='#FFFF00', zorder=0)  
        plddt_ax.axhspan(70, 90, facecolor='#6495ED', zorder=0)  
        plddt_ax.axhspan(90, 100, facecolor='#0000FF', zorder=0)   

        # Set axis labels, ticks, and limits
        plddt_ax.set(
            yticks = [0,50, 70, 90,100],
            xlim = [0,len(plddt_data)],
            ylim = [0,100],
            xlabel = 'Residue Index',
            ylabel = r'C$\regular{\alpha}$-plDDT'
        )

        return plddt_ax

    @staticmethod
    def AF3_pae_plot(pae_data:List[List]):

        pae_ax = sns.heatmap(
            pae_data, 
            cmap='Greens_r',
            vmin=0,
            vmax=32.75,
            cbar=True,
            cbar_kws={'orientation': 'horizontal', 'label':'Predicted Aligned Error (Å)', 'ticks':[0,5,10,15,20,25,30],'pad':0.2})

        pae_ax.set(
            xlabel = 'Scored Residue',
            ylabel = 'Aligned Residue',
            xticks = [0, 23, 46, 69, 92, 115, 138],
            yticks = [0, 23, 46, 69, 92, 115, 138],
            xticklabels = [1, 24, 47, 70, 93, 116, 139],
            yticklabels = [1, 24, 47, 70, 93, 116, 139],
        )

        # Rotate x-tick labels
        plt.xticks(rotation=0)

        # Remove tick marks but keep the labels
        plt.tick_params(axis='both', which='both', length=0)

        return pae_ax
