import argparse
import pandas as pd
import numpy as np
from get_distances_utils import identityEuclideanGroups

def parse_commandline():
    parser = argparse.ArgumentParser(description="Process spatial image dataset")
    
    parser.add_argument('--file', type=str, required=True, help='input csv file')
    parser.add_argument('--output', type=str, required=True, help='output csv file prefix')
    parser.add_argument('--pixelsXYZ', type=float, nargs='+', default=[0.87, 0.87, 5], help='xyz pixel dimensions')
    parser.add_argument('--cells', type=int, nargs='+', default=[0, 2, 4], help='cells to select')
    parser.add_argument('--intensity_match', type=str, nargs='+', default=['c2', 'c4', 'c0'], help='intensity match labels')

    return parser.parse_args()

def strip_whitespace(arr):
    return [x.strip() for x in arr]

def main():
    args = parse_commandline()

    input_file = args.file
    output_file = args.output
    pixelsXYZ = args.pixelsXYZ
    cells = args.cells
    intensity_match = strip_whitespace(args.intensity_match)

    dict_cells = dict(zip(cells, intensity_match))

    print(f"input = {input_file}")
    print(f"output = {output_file}")
    print(f"cells = {cells}")
    print(f"intensity_match = {intensity_match}")
    print(f"dict_cells = {dict_cells}")

    x, y, z = pixelsXYZ
    print(f"pixel dimensions x = {x}, y = {y}, z = {z}")

    df = pd.read_csv(input_file)

    for i in cells:
        ddf = df[df['cell'] == i].copy()

        df_out = identityEuclideanGroups(ddf, x, y, z, dict_cells, i)

        df_out['wmin'] = df_out['x'] - df_out['w'] / 2
        df_out['wmax'] = df_out['x'] + df_out['w'] / 2
        df_out['hmin'] = df_out['y'] - df_out['h'] / 2
        df_out['hmax'] = df_out['y'] + df_out['h'] / 2
        df_out['area'] = df_out['h'] * df_out['w']

        grouped = df_out.groupby('FinalCell').agg(
            zmin=('z', 'min'),
            zmean=('z', 'mean'),
            zmax=('z', 'max'),
            hmax=('hmax', lambda x: x.quantile(0.5)),
            hmin=('hmin', lambda x: x.quantile(0.5)),
            wmax=('wmax', lambda x: x.quantile(0.5)),
            wmin=('wmin', lambda x: x.quantile(0.5)),
            x=('x', 'mean'),
            y=('y', 'mean')
        ).reset_index()

        grouped = grouped.rename(columns={'FinalCell': 'cell_number'})
        grouped.to_csv(f"{output_file[:-4]}{i}.csv", index=False)

if __name__ == "__main__":
    main()
