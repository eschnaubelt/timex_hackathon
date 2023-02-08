from pydantic import BaseModel
from typing import (Deque, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple, Union, Type)
import ruamel.yaml


class Winding(BaseModel):
    names: List[str] = None
    currents: List[float] = None
    current_driven: List[bool] = None
    excitation_amps: List[float] = None
    heights: List[float] = None
    no_turns: List[float] = None
    inner_rs: List[float] = None
    cable_arc_lengths: List[float] = None

    sep_layer_widths: List[float] = None
    sep_layer_sigmas: List[float] = None
    sep_layer_mu_rs: List[float] = None
    sep_layer_meshed: List[bool] = None

    inn_layer_widths: List[float] = None
    inn_layer_sigmas: List[float] = None
    inn_layer_mu_rs: List[float] = None

    rs_copper_in: List[float] = None
    rs_copper_out: List[float] = None

    copper_disk: bool = None

    insulation_height: float = None

    duplicate_bnd: List[bool] = None


class Mesh(BaseModel):
    nums_elem_height_coil: List[int] = None
    nums_elem_height_cable: List[int] = None
    msh_size_from_curvature: List[int] = None
    msh_size_max: List[float] = None
    recombine_coil: List[bool] = None
    recombine_cable: List[bool] = None


class Air(BaseModel):
    height: float = None
    width: float = None


class PancakeDM(BaseModel):
    magnet_name: str = None
    windings: Winding = Winding()
    # formers: Former = Former()
    air: Air = Air()
    mesh: Mesh = Mesh()

    magnet_name: str = None         # region name


if __name__ == "__main__":

    write = True
    read = True

    def read_data(data_file_name):
        with open(data_file_name, 'r') as stream:
            yaml_str = ruamel.yaml.safe_load(stream)
        return PancakeDM(**yaml_str)

    if write:
        cct = PancakeDM()
        with open('cct_magnet_empty.yaml', 'w') as yaml_file:
            ruamel.yaml.dump(cct.dict(), yaml_file, default_flow_style=False)
    if read:
        data_file_name = 'cct1.yaml'
        data = read_data(data_file_name)
        print(data)

