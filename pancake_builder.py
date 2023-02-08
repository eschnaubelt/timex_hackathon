import gmsh

import os
import shutil
import timeit

import ruamel
import pancake_data_model as pcdm
import sys
import math
import numpy as np

def flatten(l):
    return [item for sublist in l for item in sublist]

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

class Pancake:
    def __init__(self, pancake_name, verbose=False):
        """
        Generates complete cct magnet model
        :param pancake_name: string with yaml input file name (without .yaml extension)
        :param verbose: boolean with verbose mode. True prints more to screen for debugging purpose.
        """

        self.model_name = pancake_name
        self.verbose = verbose

        self.no_subdivisions_turn = 400  # no of points per turn
        self.no_subdivisions_spline = 3  # splines are subdivided for each turn to generate nicer splines
        self.include_air = False  # include air
        self.include_cable = False  # include cable

        self.vols = None
        self.hole_vol = []
        self.coil_vols = []
        self.cable_vols = []
        self.turn_turn_lines = []
        self.thin_layers = []
        self.thin_layers_copy = []
        self.thin_layers_boundary = []
        self.thin_layers_copy_boundary = []
        self.ports = []
        self.windings_holes_bnds = []
        self.bare_bnds = []
        self.bare_hole_bnds = []
        self.windings_bnds = []
        self.outer_bnd = None
        self.aux_layer = []

        subdivide_spline = round(1 / self.no_subdivisions_spline * self.no_subdivisions_turn)

        self.firstLocalTurnTurn = 1 * subdivide_spline - 15
        self.secondLocalTurnTurn = 1 * subdivide_spline

        self.local_defect_turn = 15
        self.local_defect_idx = 2 + (self.local_defect_turn - 1) * 3

        with open(f'inputs/{self.model_name}.yaml', 'r') as stream:
            yaml_str = ruamel.yaml.safe_load(stream)
        self.pcdm = pcdm.PancakeDM(**yaml_str)

        if self.pcdm.windings.sep_layer_meshed[0]:
            self.mesh_suffix = "_explMesh"
        else:
            self.mesh_suffix = ""

        self.no_turns = self.pcdm.windings.no_turns[0]

        self.ac = 0  # action counter
        gmsh.initialize(sys.argv)

    def check_for_event(self):
        action = gmsh.onelab.getString("ONELAB/Action")
        if self.verbose:
            self.ac += 1
            # print(f"Action number {self.ac} was triggered = {action}")
        if len(action) and action[0] == "check":
            gmsh.onelab.setString("ONELAB/Action", [""])
            if self.verbose:
                print("-------------------check----------------")
            # if (gmsh.onelab.getChanged("Gmsh")):
            #     self.generate_winding()
            # if (gmsh.onelab.getChanged("GetDP")):
            #     self.load_pro()
            gmsh.fltk.update()
            gmsh.graphics.draw()
        if len(action) and action[0] == "compute":
            gmsh.onelab.setString("ONELAB/Action", [""])
            if self.verbose:
                print("-------------------compute----------------")
            #  self.generate_mesh_and_save()
            gmsh.onelab.setChanged("Gmsh", 0)
            # self.load_pro_and_solve()
            gmsh.onelab.setChanged("GetDP", 0)
            gmsh.fltk.update()
            gmsh.graphics.draw()
        return True

    def launch_interactive_gui(self):
        # if "-nopopup" not in sys.argv:
        gmsh.fltk.initialize()
        while gmsh.fltk.isAvailable() and self.check_for_event():
            gmsh.fltk.wait()
        gmsh.finalize()

    def prep(self, clear=False):
        if clear:
            if os.path.exists(self.model_name):
                shutil.rmtree(self.model_name)  # delete directory
        if not os.path.exists(self.model_name):
            os.makedirs(self.model_name)  # make new directory

    def gen_coils(self):
        def gen_coil_top_view():

            # x coordinate for the midline of the spiral
            def spiral_x_coord(a, spiral_initial_r, growth_rate):
                return (spiral_initial_r + growth_rate * a) * math.cos(a)

            # y coordinate for the midline of the spiral
            def spiral_y_coord(a, spiral_initial_r, growth_rate):
                return (spiral_initial_r + growth_rate * a) * math.sin(a)

            # differentials of spiral in x-dir
            def spiral_dx(a, spiral_initial_r, growth_rate):
                return growth_rate * math.cos(a) - (spiral_initial_r + growth_rate * a) * math.sin(a)

            # differentials of spiral in y-dir
            def spiral_dy(a, spiral_initial_r, growth_rate):
                return growth_rate * math.sin(a) + (spiral_initial_r + growth_rate * a) * math.cos(a)

            # correction to account for thickness of spiral in x-direction
            def spiral_nx(a, spiral_initial_r, growth_rate):
                return - spiral_dy(a, spiral_initial_r, growth_rate) / math.sqrt(
                    spiral_dx(a, spiral_initial_r, growth_rate) ** 2
                    + spiral_dy(a, spiral_initial_r, growth_rate) ** 2)

            # correction to account for thickness of spiral in y-direction
            def spiral_ny(a, spiral_initial_r, growth_rate):
                return spiral_dx(a, spiral_initial_r, growth_rate) / math.sqrt(
                    spiral_dx(a, spiral_initial_r, growth_rate) ** 2
                    + spiral_dy(a, spiral_initial_r, growth_rate) ** 2)

            # generate points of spiral (and cable extensions)
            def gen_spiral_points():
                # width of the cable
                w_cable = self.pcdm.windings.inn_layer_widths[winding_idx] + 2 * self.pcdm.windings.sep_layer_widths[
                    winding_idx]
                # height of the cable
                h_cable = self.pcdm.windings.heights[winding_idx]

                spiral_initial_r = self.pcdm.windings.inner_rs[
                                       winding_idx] + w_cable / 2  # inner midline radius of the spiral
                spiral_initial_theta = 0  # initial angle of spiral
                growth_rate = (self.pcdm.windings.inn_layer_widths[winding_idx] +
                               2 * self.pcdm.windings.sep_layer_widths[winding_idx]) / (
                                      2 * math.pi)  # spiral growth rate

                spiral_pnts = []  # list of all spiral points

                no_turns = self.pcdm.windings.no_turns[winding_idx]

                # loop over turns and aux regions to add cables (i.e., loop over more indices)
                for turn_idx in range(-1, math.floor(no_turns + 2)):
                    start_theta = spiral_initial_theta + turn_idx * 2 * math.pi  # start angle of this turn

                    end_theta = spiral_initial_theta + (turn_idx + 1) * 2 * math.pi  # end angle of this turn

                    # init coords for this turn
                    spiral_pnts.append([])

                    # loop of subdivision per turn and add spiral coordinates
                    for i in np.linspace(start_theta, end_theta, self.no_subdivisions_turn, endpoint=False):
                        curr_x_coord = spiral_x_coord(i, spiral_initial_r, growth_rate) + \
                                       spiral_nx(i, spiral_initial_r, growth_rate) * w_cable / 2
                        curr_y_coord = spiral_y_coord(i, spiral_initial_r, growth_rate) + \
                                       spiral_ny(i, spiral_initial_r, growth_rate) * w_cable / 2

                        spiral_pnts[turn_idx + 1].append(
                            gmsh.model.occ.addPoint(curr_x_coord, curr_y_coord, h_cable / 2))

                return spiral_pnts

            # generate inner points of spiral
            def gen_inner_points():

                inner_points_inner = []
                inner_points_outer = []

                w_ins = self.pcdm.windings.sep_layer_widths[winding_idx]
                w_cable = self.pcdm.windings.inn_layer_widths[winding_idx] + 2 * \
                          self.pcdm.windings.sep_layer_widths[
                              winding_idx]
                no_turns = self.pcdm.windings.no_turns[winding_idx]
                h_cable = self.pcdm.windings.heights[winding_idx]

                part_ins = w_ins / w_cable

                for turn_idx in range(0, math.floor(no_turns + 2)):
                    inner_points_inner.append([])
                    inner_points_outer.append([])

                    for i in range(0, self.no_subdivisions_turn):
                        inner_spiral_point = gmsh.model.getValue(0, spiral_points[turn_idx][i], [])
                        outer_spiral_point = gmsh.model.getValue(0, spiral_points[turn_idx + 1][i], [])

                        x_inner = (1 - part_ins) * inner_spiral_point[0] + part_ins * outer_spiral_point[0]
                        x_outer = part_ins * inner_spiral_point[0] + (1 - part_ins) * outer_spiral_point[0]

                        y_inner = (1 - part_ins) * inner_spiral_point[1] + part_ins * outer_spiral_point[1]
                        y_outer = part_ins * inner_spiral_point[1] + (1 - part_ins) * outer_spiral_point[1]

                        inner_points_inner[turn_idx].append(gmsh.model.occ.addPoint(x_inner, y_inner, h_cable / 2))
                        inner_points_outer[turn_idx].append(gmsh.model.occ.addPoint(x_outer, y_outer, h_cable / 2))

                return inner_points_inner, inner_points_outer

            def gen_turn_turn_lines_tsa():

                cable_arc_length = self.pcdm.windings.cable_arc_lengths[winding_idx]  # arc length of cable
                no_turns = self.pcdm.windings.no_turns[winding_idx]
                w_cable = self.pcdm.windings.inn_layer_widths[winding_idx] + 2 * self.pcdm.windings.sep_layer_widths[
                    winding_idx]
                outer_r = self.pcdm.windings.inner_rs[winding_idx] + math.floor(
                    no_turns) * w_cable  # radius of outer turn

                turn_turn_lines = []  # lines from turn to turn

                # compute indices corresponding to the cable position
                cable_in_idx = self.no_subdivisions_turn - round(cable_arc_length
                                                                 / (2 * math.pi * self.pcdm.windings.inner_rs[
                    winding_idx]) *
                                                                 self.no_subdivisions_turn)  # beginning of cable_in
                # end of last turn, i.e., beginning of cable_out
                cable_end_idx = round(no_turns % 1 * self.no_subdivisions_turn)
                # end of cable_out
                cable_out_idx = cable_end_idx + round(cable_arc_length
                                                      / (2 * math.pi * outer_r) * self.no_subdivisions_turn)

                # cable_in line
                turn_turn_lines.append(gmsh.model.occ.addLine(spiral_points[0][cable_in_idx],
                                                              spiral_points[1][cable_in_idx]))

                # turn to turn at beginning of turns
                for j in range(1, math.floor(no_turns + 2)):
                    turn_turn_lines.append(gmsh.model.occ.addLine(spiral_points[j][0], spiral_points[j + 1][0]))

                    # turn_turn_lines.append(gmsh.model.occ.addLine(spiral_points[j][self.firstLocalTurnTurn],
                    #                                               spiral_points[j + 1][self.firstLocalTurnTurn]))
                    turn_turn_lines.append(gmsh.model.occ.addLine(spiral_points[j][self.firstLocalTurnTurn],
                                                                  spiral_points[j + 1][self.firstLocalTurnTurn]))
                    turn_turn_lines.append(gmsh.model.occ.addLine(spiral_points[j][self.secondLocalTurnTurn],
                                                                  spiral_points[j + 1][self.secondLocalTurnTurn]))

                # end of last turn
                turn_turn_lines.append(
                    gmsh.model.occ.addLine(spiral_points[math.floor(no_turns + 1)][cable_end_idx],
                                           spiral_points[math.floor(no_turns + 2)][cable_end_idx]))

                # end of cable_out
                turn_turn_lines.append(
                    gmsh.model.occ.addLine(spiral_points[math.floor(no_turns + 1)][cable_out_idx],
                                           spiral_points[math.floor(no_turns + 2)][cable_out_idx]))

                return turn_turn_lines, cable_in_idx, cable_end_idx, cable_out_idx

            def gen_turn_turn_lines_meshed():

                cable_arc_length = self.pcdm.windings.cable_arc_lengths[winding_idx]  # arc length of cable
                no_turns = self.pcdm.windings.no_turns[winding_idx]
                w_cable = self.pcdm.windings.inn_layer_widths[winding_idx] + 2 * self.pcdm.windings.sep_layer_widths[
                    winding_idx]
                outer_r = self.pcdm.windings.inner_rs[winding_idx] + math.floor(
                    no_turns) * w_cable  # radius of outer turn

                turn_turn_lines = []  # lines from turn to turn

                # compute indices corresponding to the cable position
                cable_in_idx = self.no_subdivisions_turn - round(cable_arc_length
                                                                 / (2 * math.pi * self.pcdm.windings.inner_rs[
                    winding_idx]) *
                                                                 self.no_subdivisions_turn)  # beginning of cable_in
                # end of last turn, i.e., beginning of cable_out
                cable_end_idx = round(no_turns % 1 * self.no_subdivisions_turn)
                # end of cable_out
                cable_out_idx = cable_end_idx + round(cable_arc_length
                                                      / (2 * math.pi * outer_r) * self.no_subdivisions_turn)

                # cable_in line
                turn_turn_lines.append([gmsh.model.occ.addLine(inner_points_inner[0][cable_in_idx],
                                                               inner_points_outer[0][cable_in_idx]),
                                        gmsh.model.occ.addLine(inner_points_outer[0][cable_in_idx],
                                                               inner_points_inner[1][cable_in_idx])])

                # beginning of first turn
                turn_turn_lines.append([gmsh.model.occ.addLine(inner_points_inner[1][0], inner_points_outer[1][0]),
                                        gmsh.model.occ.addLine(inner_points_outer[1][0], inner_points_inner[2][0])])

                # local defect
                turn_turn_lines.append([gmsh.model.occ.addLine(inner_points_inner[self.local_defect_turn]
                                                               [self.firstLocalTurnTurn],
                                                               inner_points_outer[self.local_defect_turn]
                                                               [self.firstLocalTurnTurn]),
                                        gmsh.model.occ.addLine(inner_points_inner[self.local_defect_turn]
                                                               [self.secondLocalTurnTurn],
                                                               inner_points_outer[self.local_defect_turn]
                                                               [self.secondLocalTurnTurn])])

                # end of last turn
                turn_turn_lines.append(
                    [gmsh.model.occ.addLine(inner_points_outer[math.floor(no_turns)][cable_end_idx],
                                            inner_points_inner[math.floor(no_turns + 1)][cable_end_idx]),
                     gmsh.model.occ.addLine(inner_points_inner[math.floor(no_turns + 1)][cable_end_idx],
                                            inner_points_outer[math.floor(no_turns + 1)][cable_end_idx])])

                # end of cable_out
                turn_turn_lines.append(
                    [gmsh.model.occ.addLine(inner_points_outer[math.floor(no_turns)][cable_out_idx],
                                            inner_points_inner[math.floor(no_turns + 1)][cable_out_idx]),
                     gmsh.model.occ.addLine(inner_points_inner[math.floor(no_turns + 1)][cable_out_idx],
                                            inner_points_outer[math.floor(no_turns + 1)][cable_out_idx])])

                return turn_turn_lines, cable_in_idx, cable_end_idx, cable_out_idx

            def gen_splines_tsa():

                # fraction of incomplete turn
                incomplete_turn = round((self.pcdm.windings.no_turns[winding_idx] % 1) * self.no_subdivisions_turn)

                subdivide_spline = round(1 / self.no_subdivisions_spline * self.no_subdivisions_turn)

                # collect all indices of points that need to be connected by splintes
                connecting_idx = [i * subdivide_spline for i in range(0, self.no_subdivisions_spline)] + \
                                 [incomplete_turn, cable_in_index, cable_out_index] + \
                                 [self.firstLocalTurnTurn]

                connecting_idx = list(set(connecting_idx))
                connecting_idx.sort()

                # cable_in inner spline
                splines = [[gmsh.model.occ.addSpline(spiral_points[0][cable_in_index:] +
                                                     [spiral_points[1][0]])]]

                # inner turn splines (only when thin shell is used)
                for j in range(1, math.floor(self.pcdm.windings.no_turns[winding_idx] + 2)):
                    splines.append([])

                    # splines inside turn
                    for i in range(0, len(connecting_idx) - 1):
                        splines[j].append(gmsh.model.occ.addSpline(spiral_points[j]
                                                                   [connecting_idx[i]:connecting_idx[i + 1] + 1]))

                    # last spline connects to next turn
                    splines[j].append(gmsh.model.occ.addSpline(spiral_points[j]
                                                               [connecting_idx[i + 1]:] + [spiral_points[j + 1]
                                                                                           [0]]))

                splines.append([])

                # spline for last (possibly incomplete) turn
                no_spline_before_cable_out = self.no_subdivisions_turn

                for i in range(0, len(connecting_idx) - 1):
                    if connecting_idx[i] > cable_end_index:
                        no_spline_before_cable_out = i - 1
                        break
                    splines[-1].append(gmsh.model.occ.addSpline(spiral_points[-1]
                                                                [connecting_idx[i]:connecting_idx[i + 1] + 1]))

                self.splines = splines

                return splines, no_spline_before_cable_out, connecting_idx

            def gen_splines_meshed():

                # fraction of incomplete turn
                incomplete_turn = round((self.pcdm.windings.no_turns[winding_idx] % 1) * self.no_subdivisions_turn)

                subdivide_spline = round(1 / self.no_subdivisions_spline * self.no_subdivisions_turn)

                # collect all indices of points that need to be connected by splintes
                connecting_idx = [i * subdivide_spline for i in range(0, self.no_subdivisions_spline)] + \
                                 [incomplete_turn, cable_in_index, cable_out_index] + \
                                 [self.firstLocalTurnTurn]

                connecting_idx = list(set(connecting_idx))
                connecting_idx.sort()

                # cable_in inner spline
                splines_inner = [[gmsh.model.occ.addSpline(inner_points_inner[0][cable_in_index:] +
                                                           [inner_points_inner[1][0]])]]
                splines_outer = [[gmsh.model.occ.addSpline(inner_points_outer[0][cable_in_index:] +
                                                           [inner_points_outer[1][0]])]]

                # inner turn splines (only when thin shell is used)
                for j in range(1, math.floor(self.pcdm.windings.no_turns[winding_idx] + 1)):
                    splines_inner.append([])
                    splines_outer.append([])

                    # splines inside turn
                    for i in range(0, len(connecting_idx) - 1):
                        splines_inner[j].append(gmsh.model.occ.addSpline(inner_points_inner[j]
                                                                         [connecting_idx[i]:connecting_idx[i + 1] + 1]))

                        splines_outer[j].append(gmsh.model.occ.addSpline(inner_points_outer[j]
                                                                         [connecting_idx[i]:connecting_idx[i + 1] + 1]))

                    splines_inner[j].append(gmsh.model.occ.addSpline(inner_points_inner[j]
                                                                     [connecting_idx[i + 1]:] + [
                                                                         inner_points_inner[j + 1]
                                                                         [0]]))

                    splines_outer[j].append(gmsh.model.occ.addSpline(inner_points_outer[j]
                                                                     [connecting_idx[i + 1]:] + [
                                                                         inner_points_outer[j + 1]
                                                                         [0]]))

                splines_inner.append([])
                splines_outer.append([])

                # spline for last (possibly incomplete) turn
                no_spline_before_cable_out = self.no_subdivisions_turn

                for i in range(0, len(connecting_idx) - 1):
                    if connecting_idx[i] > cable_end_index:
                        no_spline_before_cable_out = i - 1
                        break
                    splines_inner[-1].append(gmsh.model.occ.addSpline(inner_points_inner[-1][connecting_idx[i]:
                                                                                             connecting_idx[
                                                                                                 i + 1] + 1]))
                    splines_outer[-1].append(gmsh.model.occ.addSpline(inner_points_outer[-1][connecting_idx[i]:
                                                                                             connecting_idx[
                                                                                                 i + 1] + 1]))

                return splines_inner, splines_outer, no_spline_before_cable_out, connecting_idx

            def gen_surfaces_tsa():

                # cable_in surface
                wire_tags = [gmsh.model.occ.addWire([turn_turn_lines[0]] + splines[0] +
                                                    splines[1][:-1]),
                             gmsh.model.occ.addWire([turn_turn_lines[0]] + splines[0] +
                                                    [turn_turn_lines[1], splines[1][-1]])]

                surfaces = [gmsh.model.occ.addPlaneSurface([wire_tags[0]]),
                            gmsh.model.occ.addPlaneSurface([wire_tags[1]])]

                for j in range(1, math.floor(self.pcdm.windings.no_turns[winding_idx] + 2)):
                    for i in range(0, 2):
                        wire_tags.append(gmsh.model.occ.addWire([turn_turn_lines[j + i + (j - 1) * 2], splines[j][i],
                                                                 turn_turn_lines[j + i + (j - 1) * 2 + 1],
                                                                 splines[j + 1][i]]))

                        surfaces.append(gmsh.model.occ.addPlaneSurface([wire_tags[-1]]))

                    i = 2

                    # last turn or not (last one possibly incomplete)
                    if j < math.floor(self.pcdm.windings.no_turns[winding_idx] + 1):
                        wire_tags.append(
                            gmsh.model.occ.addWire([turn_turn_lines[j + i + (j - 1) * 2]] + splines[j][i::] +
                                                   [turn_turn_lines[j + i + (j - 1) * 2 + 1]] +
                                                   splines[j + 1][i::]))

                    else:
                        wire_tags.append(gmsh.model.occ.addWire(
                            [turn_turn_lines[j + i + (j - 1) * 2]] + splines[j][i::no_spline_before_cable_out] +
                            [turn_turn_lines[j + i + (j - 1) * 2 + 1]] +
                            splines[j + 1][i::no_spline_before_cable_out]))

                    surfaces.append(gmsh.model.occ.addPlaneSurface([wire_tags[-1]]))

                # cable_out
                wire_tags.append(gmsh.model.occ.addWire([turn_turn_lines[-2],
                                                         splines[-2][no_spline_before_cable_out],
                                                         turn_turn_lines[-1],
                                                         splines[-1][no_spline_before_cable_out]]))

                surfaces.append(gmsh.model.occ.addPlaneSurface([wire_tags[-1]]))

                return surfaces

            def gen_surfaces_meshed():

                # hole, cable in surfaces
                wire_tags = [
                    gmsh.model.occ.addWire(turn_turn_lines[0] + [splines_inner[0][0]] + splines_inner[1][0:-1]),
                    gmsh.model.occ.addWire([turn_turn_lines[0][0], splines_inner[0][0],
                                            turn_turn_lines[1][0], splines_outer[0][0]])]

                surfaces = []

                # FIXME HERE!
                # inner surfaces
                first_wire = [spline for list_of_splines in splines_inner[1:self.local_defect_turn] for
                              spline in list_of_splines] + \
                             [splines_inner[self.local_defect_turn][0]] + \
                             [turn_turn_lines[1][0]] + \
                             [spline for list_of_splines in splines_outer[1:self.local_defect_turn] for
                              spline in list_of_splines] + \
                             [splines_outer[self.local_defect_turn][0]] + \
                             [turn_turn_lines[2][0]]

                wire_tags.append(gmsh.model.occ.addWire(first_wire))

                second_wire = [splines_inner[self.local_defect_turn][1],
                               turn_turn_lines[2][0],
                               splines_outer[self.local_defect_turn][1],
                               turn_turn_lines[2][1]]

                wire_tags.append(gmsh.model.occ.addWire(second_wire))

                third_wire = splines_inner[self.local_defect_turn][2:] + \
                             [spline for list_of_splines in splines_inner[self.local_defect_turn + 1:-1] for
                              spline in list_of_splines] + \
                             splines_inner[-1][:no_spline_before_cable_out] + \
                             [turn_turn_lines[2][1]] + \
                             splines_outer[self.local_defect_turn][2:] + \
                             [spline for list_of_splines in splines_outer[self.local_defect_turn + 1:-1] for
                              spline in list_of_splines] + \
                             splines_outer[-1][:no_spline_before_cable_out] + \
                             [turn_turn_lines[3][1]]

                wire_tags.append(gmsh.model.occ.addWire(third_wire))

                # insulation surface
                unpack_outer_splines = [spline for list_of_splines in splines_outer[:-2] for spline in list_of_splines] \
                                       + splines_outer[-2][:no_spline_before_cable_out + 1]

                unpack_inner_splines = [spline for list_of_splines in splines_inner[2:] for spline in list_of_splines]

                wire_tags.append(gmsh.model.occ.addWire([turn_turn_lines[0][1], splines_inner[1][-1]] +
                                                        unpack_inner_splines + [turn_turn_lines[-1][0]] +
                                                        unpack_outer_splines))

                # cable out surface
                wire_tags.append(gmsh.model.occ.addWire([turn_turn_lines[-1][1], splines_inner[-1][-1],
                                                         turn_turn_lines[-2][1], splines_outer[-1][-1]]))

                for wire in wire_tags[:]:
                    surfaces.append(gmsh.model.occ.addPlaneSurface([wire]))

                return surfaces

            def remove_points():
                for point_list in spiral_points:
                    for point in point_list:
                        gmsh.model.occ.remove([(0, point)])

            def remove_inner_points():
                for point_list in inner_points_outer:
                    for point in point_list:
                        gmsh.model.occ.remove([(0, point)])

                for point_list in inner_points_inner:
                    for point in point_list:
                        gmsh.model.occ.remove([(0, point)])

            # generate closely packed points on spirals with which splines will be defined
            spiral_points = gen_spiral_points()

            if self.pcdm.windings.sep_layer_meshed[winding_idx]:  # explicitly mesh separation layer
                gmsh.model.occ.synchronize()
                inner_points_inner, inner_points_outer = gen_inner_points()
                # generate lines between turns
                turn_turn_lines, cable_in_index, cable_end_index, cable_out_index = gen_turn_turn_lines_meshed()
                # generate splines for the turn-to-turn interfaces
                splines_inner, splines_outer, no_spline_before_cable_out, connecting_idx = gen_splines_meshed()
                surfaces = gen_surfaces_meshed()
                remove_inner_points()
            else:  # thin shell approach
                # generate lines between turns
                turn_turn_lines, cable_in_index, cable_end_index, cable_out_index = gen_turn_turn_lines_tsa()
                # generate splines for the turn-to-turn interfaces
                splines, no_spline_before_cable_out, connecting_idx = gen_splines_tsa()
                surfaces = gen_surfaces_tsa()

            remove_points()

            return surfaces, turn_turn_lines

        def gen_coil():
            def extrude_surfaces():

                if self.include_air:
                    extr_hole = gmsh.model.occ.extrude([(2, surfaces[0])], 0, 0,
                                                       -self.pcdm.windings.heights[winding_idx])
                else:
                    extr_hole = []

                    gmsh.model.occ.remove([(2, surfaces[0])])

                surface_dim_tags = [(2, surface) for surface in surfaces[1:]]
                extr_coil = gmsh.model.occ.extrude(surface_dim_tags, 0, 0, -self.pcdm.windings.heights[winding_idx],
                                                   numElements=[self.pcdm.mesh.nums_elem_height_coil[winding_idx]],
                                                   recombine=self.pcdm.mesh.recombine_coil[winding_idx])

                cable_dim_tags = [(2, surfaces[1]), (2, surfaces[-1])]
                h_cable = self.pcdm.air.height / 2 - self.pcdm.windings.heights[winding_idx] / 2
                extr_cable = gmsh.model.occ.extrude(cable_dim_tags, 0, 0, h_cable,
                                                    numElements=[self.pcdm.mesh.nums_elem_height_cable[winding_idx]],
                                                    recombine=self.pcdm.mesh.recombine_cable[winding_idx])

                self.ports.append([extr_cable[0][1], extr_cable[6][1]])

                return extr_hole, extr_coil, extr_cable

            extruded_hole, extruded_coil, extruded_cable = extrude_surfaces()
            return extruded_hole, extruded_coil, extruded_cable

        for winding_idx in range(len(self.pcdm.windings.names)):
            surfaces, turn_turn_lines = gen_coil_top_view()

            extruded_hole, extruded_coil, extruded_cable = gen_coil()

            self.turn_turn_lines.append(turn_turn_lines)
            self.cable_vols.append([(dim, tag) for (dim, tag) in extruded_cable if dim == 3])
            self.coil_vols.append([(dim, tag) for (dim, tag) in extruded_coil if dim == 3])
            self.hole_vol.append([(dim, tag) for (dim, tag) in extruded_hole if dim == 3])

            if self.pcdm.windings.sep_layer_meshed[winding_idx]:
                self.defect_vols = [self.coil_vols[0].pop(2)]
            else:

                self.coilSort = [[], [], [], []]

                self.localDefectSort = [[], [], [], []]

                self.coilSort[0] = self.coil_vols[winding_idx][::2 * self.no_subdivisions_spline]
                self.coilSort[1] = self.coil_vols[winding_idx][1::2 * self.no_subdivisions_spline] + \
                                   self.coil_vols[winding_idx][2::2 * self.no_subdivisions_spline]
                self.coilSort[2] = self.coil_vols[winding_idx][3::2 * self.no_subdivisions_spline]
                self.coilSort[3] = self.coil_vols[winding_idx][4::2 * self.no_subdivisions_spline] +\
                                   self.coil_vols[winding_idx][5::2 * self.no_subdivisions_spline]

                self.defect_vols = [self.coil_vols[winding_idx].pop(self.local_defect_idx)]

                for i in range(len(self.coilSort)):
                    try:
                        self.coilSort[i].remove(self.defect_vols[0])
                        self.localDefectSort[i].append(self.defect_vols[0])
                    except ValueError:
                        pass  # do nothing!

    def gen_air(self):
        h_air = self.pcdm.air.height
        w_air = self.pcdm.air.width

        # collect all holes, coils and cables
        objects = [dimTag for list_of_dimTags in self.hole_vol + self.coil_vols + self.cable_vols for dimTag in
                   list_of_dimTags]

        if self.include_air:
            # outer air box
            outer_box = gmsh.model.occ.addBox(-w_air, -w_air, -h_air / 2, 2 * w_air, 2 * w_air, h_air)

            self.vols, _ = gmsh.model.occ.fragment([(3, outer_box)], objects)
        else:
            # self.vols, _ = gmsh.model.occ.fragment([objects[0]], objects[1:])
            self.vols = objects

        gmsh.model.occ.synchronize()

    def get_and_copy_thin_layers(self):
        for winding_idx in range(len(self.pcdm.windings.names)):

            if not self.pcdm.windings.sep_layer_meshed[winding_idx]:

                coil_boundary_combined = gmsh.model.getBoundary(self.coil_vols[winding_idx] +
                                                                self.defect_vols, oriented=False)
                coil_boundary_uncombined = gmsh.model.getBoundary(self.coil_vols[winding_idx] +
                                                                  self.defect_vols, combined=False,
                                                                  oriented=False)

                hole_boundary = gmsh.model.getBoundary(self.hole_vol[winding_idx], oriented=False)

                turn_turn_surfaces = []

                for line in self.turn_turn_lines[winding_idx]:
                    turn_turn_surfaces.extend(gmsh.model.get_adjacencies(1, line)[0])

                # remove duplicate entries in thin layer while preserving order

                thin_layer_list = [(dim, tag) for (dim, tag) in coil_boundary_uncombined if (dim, tag) not in
                                   coil_boundary_combined and tag not in turn_turn_surfaces]

                thin_layer_list.sort()

                thin_layer_list = unique(thin_layer_list)

                self.thin_layers.append(thin_layer_list)

                thin_layer_copy = [] # gmsh.model.occ.copy(self.thin_layers[winding_idx])

                self.thin_layers_copy.append(thin_layer_copy)

                # self.thin_layers_copy.append(gmsh.model.occ.fragment([thin_layer_copy[0]], thin_layer_copy[1:])[0])

                gmsh.model.occ.synchronize()

                self.thin_layers_boundary.append(gmsh.model.getBoundary(self.thin_layers[winding_idx], oriented=False))
                self.thin_layers_copy_boundary.append(gmsh.model.getBoundary(self.thin_layers_copy[winding_idx],
                                                                             oriented=False))

                # # same mesh for copy of thin layer (not automatically the case)
                # gmsh.model.mesh.setPeriodic(2, [tag for (dim, tag) in self.thin_layers_copy[winding_idx]],
                #                             [tag for (dim, tag) in self.thin_layers[winding_idx]],
                #                             [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])


                # FIXME: GENERALIZE
                self.layerSort = [[], [], [], []]
                self.layerSort[0] = [self.thin_layers[winding_idx][0]] + flatten([self.thin_layers[winding_idx]
                                                                          [10 + 14 * i:15 + 14 * i] for i in
                                                                          range(math.ceil(self.no_turns) + 1)])

                self.layerSort[1] = self.thin_layers[winding_idx][1::14] + self.thin_layers[winding_idx][2::14]
                self.layerSort[2] = flatten([self.thin_layers[winding_idx][3 + 14 * i:8 + 14 * i] for i in
                                                                          range(math.ceil(self.no_turns) + 1)])

                self.layerSort[3] = self.thin_layers[winding_idx][8::14] + self.thin_layers[winding_idx][9::14]

                # FIXME: this seems to be incorrect
                self.aux_layer.append(
                    list(set([(dim, tag) for (dim, tag) in coil_boundary_combined if tag not in turn_turn_surfaces])))

    def gen_physical_groups(self):
        def get_coil_bnd():
            winding_vols = [dimTag for list_of_dimTags in self.cable_vols + self.coil_vols for dimTag in
                            list_of_dimTags]

            ports_unrolled = [ports for list_of_ports in self.ports for ports in list_of_ports]

            self.windings_bnds = [(dim, tag) for (dim, tag) in gmsh.model.getBoundary(winding_vols, oriented=False) if
                                  tag not in ports_unrolled]

        def get_coil_hole_bnd():
            winding_hole_vols = [dimTag for list_of_dimTags in self.cable_vols + self.coil_vols + self.hole_vol
                                 for dimTag in list_of_dimTags]

            ports_unrolled = [ports for list_of_ports in self.ports for ports in list_of_ports]

            self.windings_holes_bnds = [(dim, tag) for (dim, tag) in gmsh.model.getBoundary(winding_hole_vols,
                                                                                            oriented=False)
                                        if tag not in ports_unrolled]

        def add_bare_bnd():
            bare_cable_vols = [dimTag for list_of_dimTags in self.cable_vols for dimTag in list_of_dimTags] + \
                              [self.coil_vols[winding_idx][0], self.coil_vols[winding_idx][1],
                               self.coil_vols[winding_idx][-1]]

            ports_unrolled = [ports for list_of_ports in self.ports for ports in list_of_ports]

            self.bare_bnds.extend([(dim, tag) for (dim, tag) in gmsh.model.getBoundary(bare_cable_vols,
                                                                                       oriented=False)
                                   if tag not in ports_unrolled])

        def add_bare_hole_bnd():

            bare_cable_hole_vols = [dimTag for list_of_dimTags in self.cable_vols + self.hole_vol for dimTag
                                    in list_of_dimTags] + \
                                   [self.coil_vols[winding_idx][0], self.coil_vols[winding_idx][1],
                                    self.coil_vols[winding_idx][-1]]

            ports_unrolled = [ports for list_of_ports in self.ports for ports in list_of_ports]

            self.bare_hole_bnds = [(dim, tag) for (dim, tag) in gmsh.model.getBoundary(bare_cable_hole_vols,
                                                                                       oriented=False)
                                   if tag not in ports_unrolled]

        def get_outer_bnd():
            self.outer_bnd = [(dim, tag) for (dim, tag) in gmsh.model.getBoundary([self.vols[-1]], oriented=False)
                              if (dim, tag) not in self.windings_holes_bnds and tag not in self.ports]

        def add_physical_groups_windings_tsa():

            for i in range(len(self.coilSort)):
                gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.coilSort[i]], tag=i+20)
                gmsh.model.setPhysicalName(dim=3, tag=i+20, name="Coil " + str(i))

                gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.localDefectSort[i]], tag=i+40)
                gmsh.model.setPhysicalName(dim=3, tag=i+40, name="Local Def " + str(i))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.defect_vols], tag=999)
            gmsh.model.setPhysicalName(dim=3, tag=999, name="Local Defect  " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.cable_vols[winding_idx]], tag=2)
            gmsh.model.setPhysicalName(dim=3, tag=2, name="Cable " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.hole_vol[winding_idx]], tag=4)
            gmsh.model.setPhysicalName(dim=3, tag=4, name="Hole " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.thin_layers_copy[winding_idx]],
                                        tag=13)
            gmsh.model.setPhysicalName(dim=2, tag=13, name="Layer Copy " + str(winding_idx))

            for i in range(len(self.layerSort)):
                gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.layerSort[i]], tag=i + 2000)
                gmsh.model.setPhysicalName(dim=2, tag=i + 2000, name="Layer " + str(i))

            gmsh.model.addPhysicalGroup(dim=2, tags=[self.ports[winding_idx][0]], tag=20)
            gmsh.model.setPhysicalName(dim=2, tag=20, name="Port 0 " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[self.ports[winding_idx][1]], tag=21)
            gmsh.model.setPhysicalName(dim=2, tag=21, name="Port 1 " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2,
                                        tags=[tag for (dim, tag) in self.aux_layer[winding_idx]],
                                        tag=1000)
            gmsh.model.setPhysicalName(dim=2, tag=1000, name="Auxiliary Group" + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.thin_layers[winding_idx]], tag=10)
            gmsh.model.setPhysicalName(dim=2, tag=10, name="Layer Res Before Crack " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=1, tags=[tag for (dim, tag) in self.thin_layers_boundary[winding_idx]],
                                        tag=100)
            gmsh.model.setPhysicalName(dim=1, tag=100, name="Layer Boundary " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=1,
                                        tags=[tag for (dim, tag) in self.thin_layers_copy_boundary[winding_idx]],
                                        tag=101)
            gmsh.model.setPhysicalName(dim=1, tag=101, name="Layer Copy Boundary " + str(winding_idx))

        def add_physical_groups_windings_meshed():
            gmsh.model.addPhysicalGroup(dim=3,
                                        tags=[self.coil_vols[winding_idx][1][1], self.coil_vols[winding_idx][0][1],
                                              self.coil_vols[winding_idx][2][1],
                                              self.coil_vols[winding_idx][-1][1]], tag=1)
            gmsh.model.setPhysicalName(dim=3, tag=1, name="Bare " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[self.coil_vols[winding_idx][3][1]], tag=5)
            gmsh.model.setPhysicalName(dim=3, tag=5, name="Insulation " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.defect_vols], tag=999)
            gmsh.model.setPhysicalName(dim=3, tag=999, name="Local Defect  " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.cable_vols[winding_idx]], tag=2)
            gmsh.model.setPhysicalName(dim=3, tag=2, name="Cable " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=3, tags=[tag for (dim, tag) in self.hole_vol[winding_idx]], tag=4)
            gmsh.model.setPhysicalName(dim=3, tag=4, name="Hole " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[self.ports[winding_idx][0]], tag=20)
            gmsh.model.setPhysicalName(dim=2, tag=20, name="Port 0 " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[self.ports[winding_idx][1]], tag=21)
            gmsh.model.setPhysicalName(dim=2, tag=21, name="Port 1 " + str(winding_idx))

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.bare_bnds], tag=19)
            gmsh.model.setPhysicalName(dim=2, tag=19, name="Bare Cable Bnd")

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.bare_hole_bnds], tag=16)
            gmsh.model.setPhysicalName(dim=2, tag=16, name="Bare Hole Cable Bnd")

        def add_physical_groups():
            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.outer_bnd], tag=30)
            gmsh.model.setPhysicalName(dim=2, tag=30, name="outer_bnd")

            gmsh.model.addPhysicalGroup(dim=3, tags=[self.vols[-1][1]], tag=3)
            gmsh.model.setPhysicalName(dim=3, tag=3, name="Air")

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.windings_bnds], tag=17)
            gmsh.model.setPhysicalName(dim=2, tag=17, name="Coil Cable Bnd")

            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.windings_holes_bnds], tag=18)
            gmsh.model.setPhysicalName(dim=2, tag=18, name="Coil Cable Hole Bnd")

        get_coil_bnd()
        get_coil_hole_bnd()
        get_outer_bnd()

        for winding_idx in range(len(self.pcdm.windings.names)):

            if self.pcdm.windings.sep_layer_meshed[winding_idx]:
                add_bare_bnd()
                add_bare_hole_bnd()
                add_physical_groups_windings_meshed()
            else:
                add_physical_groups_windings_tsa()

            add_physical_groups()

            if not self.include_air:
                gmsh.model.occ.remove(self.hole_vol[winding_idx], recursive=True)

                self.hole_vol = [()]

            if not self.include_cable:
                gmsh.model.occ.remove(self.cable_vols[winding_idx], recursive=True)

                self.cable_vols = [()]

    def gen_full_geom(self):
        if self.verbose:
            print(f'Preparing folders')
        self.prep(clear=False)

        # self.assemble_pro()

        start_time = timeit.default_timer()
        if self.verbose:
            print(f'Started generating geometry')

        self.gen_coils()
        self.gen_air()
        self.get_and_copy_thin_layers()

        self.gen_physical_groups()
        # self.gen_regions_file()
        gmsh.model.occ.synchronize()
        if self.verbose:
            print(f'Generating geometry took {timeit.default_timer() - start_time:.2f} s')

    def generate_mesh_without_crack(self):
        # for winding_idx in range(len(self.pcdm.windings.names)):
        #     if not self.pcdm.windings.sep_layer_meshed[winding_idx]:
        #         print(self.splines)
        #         for turn_idx in range(len(self.splines)):
        #             mesh_length = 0.001
        #             w_cable = self.pcdm.windings.inn_layer_widths[winding_idx] + 2 * self.pcdm.windings.sep_layer_widths[
        #                     winding_idx]
        #             big_spline_transMesh = round((self.pcdm.windings.inner_rs[winding_idx] + turn_idx * w_cable)/
        #                                          mesh_length)
        #
        #             spline_per_turn = self.splines[turn_idx]
        #             for i in range(len(spline_per_turn)):
        #                 if i == 0 or i == 2 or i == 4 or i == 5 or i == 6:
        #                     gmsh.model.mesh.setTransfiniteCurve(spline_per_turn[i], big_spline_transMesh)
        #                     pass
        #                 elif i == 1 or i == 3:
        #                     gmsh.model.mesh.setTransfiniteCurve(spline_per_turn[i], 2 + turn_idx)
        #
        #                 else:
        #                     gmsh.model.mesh.setTransfiniteCurve(spline_per_turn[i], 2 + turn_idx)


        # self.s

        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        # math_eval_field = gmsh.model.mesh.field.add("MathEval")
        # gmsh.model.mesh.field.setString(math_eval_field, "F", "0.0000001/(x^2 + y^2)^5")

        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", self.pcdm.mesh.msh_size_from_curvature[0])
        # gmsh.option.setNumber("Mesh.MeshSizeFromCurvatureIsotropic", self.pcdm.mesh.msh_size_from_curvature[0])

        # FIXME: maybe use Mesh.MeshSizeFactor instead
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.00000001)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.pcdm.mesh.msh_size_max[0])
        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        # gmsh.option.setNumber("Mesh.Algorithm3D", 0)
        # gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.option.setNumber("Mesh.Optimize", 1)
        # gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)

        gmsh.option.setNumber("Mesh.QualityType", 2)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        # gmsh.option.setNumber("Mesh.AngleToleranceFacetOverlap", 0.2)

        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 0)

        # Mesh.RecombinationAlgorithm
        #
        # for tag in [tag for (dim, tag) in gmsh.model.occ.getEntities(dim=2)]:
        #     gmsh.model.mesh.setTransfiniteSurface(tag)

        try:
          # self.launch_interactive_gui()
          gmsh.model.mesh.generate()
        except:
          gmsh.finalize()
          raise Exception("Could not create mesh!")

        # _, eleTags, _ = gmsh.model.mesh.getElements(dim=3)

        # q_sicn = gmsh.model.mesh.getElementQualities(eleTags[0], "minSICN")
        # q_sige = gmsh.model.mesh.getElementQualities(eleTags[0], "minSIGE")
        # q_gamma = gmsh.model.mesh.getElementQualities(eleTags[0], "gamma")
        #
        # with open(self.model_name + "/mesh_quality" + self.mesh_suffix + ".txt", 'w+', encoding='utf-8') as f:
        #     f.write("sicn_min = " + str(min(q_sicn)) + "\n")
        #     f.write("sige_min = " + str(min(q_sige)) + "\n")
        #     f.write("gamma_min = " + str(min(q_gamma)) + "\n")
        #
        #     f.write("sicn_max = " + str(max(q_sicn)) + "\n")
        #     f.write("sige_max = " + str(max(q_sige)) + "\n")
        #     f.write("gamma_max = " + str(max(q_gamma)) + "\n")
        #
        #     f.write("sicn_avg = " + str(sum(q_sicn) / len(q_sicn)) + "\n")
        #     f.write("sige_avg = " + str(sum(q_sige) / len(q_sige)) + "\n")
        #     f.write("gamma_avg = " + str(sum(q_gamma) / len(q_gamma)) + "\n")

        gmsh.write(self.model_name + "/" + self.model_name + self.mesh_suffix + "_noCrack.msh")

    def use_crack_plugin(self):
        if not self.pcdm.windings.sep_layer_meshed[0]:
            gmsh.plugin.setNumber("Crack", "Dimension", 2)
            gmsh.plugin.setNumber("Crack", "PhysicalGroup", 10)

            gmsh.plugin.setNumber("Crack", "OpenBoundaryPhysicalGroup", 100)
            gmsh.plugin.setNumber("Crack", "NewPhysicalGroup", 12)
            gmsh.plugin.setNumber("Crack", "NormalZ", -1)
            gmsh.plugin.setNumber("Crack", "DebugView", 1)
            gmsh.plugin.setNumber("Crack", "AuxiliaryPhysicalGroup", 0)

            gmsh.plugin.run("Crack")

            # FIXME: General no of windings
            gmsh.model.addPhysicalGroup(dim=2, tags=[tag for (dim, tag) in self.thin_layers[0]], tag=11)
            gmsh.model.setPhysicalName(dim=2, tag=11, name="Layer Down After Crack")
            gmsh.model.setPhysicalName(dim=2, tag=12, name="Layer Up After Crack")

            # WORKAROUND UNTIL CRACK PLUGIN FIXED
            # i = 1000
            # layer_down_after_crack = []

            # open_bnd = gmsh.model.getBoundary(self.thin_layers[0], oriented=False, combined=False)

            # gmsh.model.addPhysicalGroup(dim=1, tags=[tag for (dim, tag) in open_bnd],
            # tag=200)
            # gmsh.model.setPhysicalName(dim=1, tag=200, name="Layer Boundary Open")

            # for sur in self.thin_layers[0]:
            # gmsh.model.addPhysicalGroup(dim=2, tags=[sur[1]], tag=i)
            # gmsh.model.setPhysicalName(dim=2, tag=i, name="Layer Res Before Crack" + str(i))

            # gmsh.plugin.setNumber("Crack", "Dimension", 2)
            # gmsh.plugin.setNumber("Crack", "PhysicalGroup", i)
            # gmsh.plugin.setNumber("Crack", "OpenBoundaryPhysicalGroup", 200)
            # gmsh.plugin.setNumber("Crack", "NewPhysicalGroup", i + 1000)
            # gmsh.plugin.run("Crack")

            # gmsh.model.setPhysicalName(dim=2, tag=i, name="Layer Up After Crack" + str(i))

            # layer_down_after_crack.append(i)

            # i = i + 1

            gmsh.write(self.model_name + "/" + self.model_name + self.mesh_suffix + ".msh")

    def generate_cuts(self):

        gmsh.plugin.setString("HomologyComputation", "DomainPhysicalGroups", "3")
        gmsh.plugin.setString("HomologyComputation", "SubdomainPhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "ReductionImmunePhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "DimensionOfChainsToSave", "1")
        gmsh.plugin.setString("HomologyComputation", "Filename", "Homology")
        gmsh.plugin.setNumber("HomologyComputation", "ComputeHomology", 0)
        gmsh.plugin.setNumber("HomologyComputation", "ComputeCohomology", 1)
        gmsh.plugin.setNumber("HomologyComputation", "CreatePostProcessingViews", 1)
        gmsh.plugin.run("HomologyComputation")

        gmsh.plugin.setString("HomologyComputation", "DomainPhysicalGroups", "3, 4")
        gmsh.plugin.setString("HomologyComputation", "SubdomainPhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "ReductionImmunePhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "DimensionOfChainsToSave", "1")
        gmsh.plugin.setString("HomologyComputation", "Filename", "Homology")
        gmsh.plugin.setNumber("HomologyComputation", "ComputeHomology", 0)
        gmsh.plugin.setNumber("HomologyComputation", "ComputeCohomology", 1)
        gmsh.plugin.setNumber("HomologyComputation", "CreatePostProcessingViews", 1)
        gmsh.plugin.run("HomologyComputation")

        gmsh.plugin.setString("HomologyComputation", "DomainPhysicalGroups", "3, 4, 5")
        gmsh.plugin.setString("HomologyComputation", "SubdomainPhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "ReductionImmunePhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "DimensionOfChainsToSave", "1")
        gmsh.plugin.setString("HomologyComputation", "Filename", "Homology")
        gmsh.plugin.setNumber("HomologyComputation", "ComputeHomology", 0)
        gmsh.plugin.setNumber("HomologyComputation", "ComputeCohomology", 1)
        gmsh.plugin.setNumber("HomologyComputation", "CreatePostProcessingViews", 1)
        gmsh.plugin.run("HomologyComputation")

        gmsh.plugin.setString("HomologyComputation", "DomainPhysicalGroups", "3, 5")
        gmsh.plugin.setString("HomologyComputation", "SubdomainPhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "ReductionImmunePhysicalGroups", "")
        gmsh.plugin.setString("HomologyComputation", "DimensionOfChainsToSave", "1")
        gmsh.plugin.setString("HomologyComputation", "Filename", "Homology")
        gmsh.plugin.setNumber("HomologyComputation", "ComputeHomology", 0)
        gmsh.plugin.setNumber("HomologyComputation", "ComputeCohomology", 1)
        gmsh.plugin.setNumber("HomologyComputation", "CreatePostProcessingViews", 1)
        gmsh.plugin.run("HomologyComputation")

        gmsh.write(self.model_name + "/" + self.model_name + self.mesh_suffix + ".msh")

    def create_model(self):
        self.gen_full_geom()
        self.generate_mesh_without_crack()
        self.use_crack_plugin()
        gmsh.finalize()


if __name__ == "__main__":
    model_name = 'thermal_coil_TSA'
    with_GUI = True
    verbose = True

    pancake = Pancake(model_name, verbose=verbose)
    # cct.prep(clear=True)
    # cct.generate_pro()

    # cct.gen_full_geom()
    # cct.generate_cuts()
    # cct.generate_mesh_and_save()
    # cct.assemble_pro()
    # cct.check_solve()

    if with_GUI:
        # pancake.prep(clear=True)
        pancake.gen_full_geom()
        pancake.generate_mesh_without_crack()
        pancake.use_crack_plugin()
        # pancake.generate_cuts()
    # pancake.assemble_pro()
    # pancake.load_pro()
    # pancake.check_solve()
    else:
        pancake.gen_full_geom()
        pancake.generate_mesh_without_crack()
        pancake.use_crack_plugin()

    if with_GUI:

        pancake.launch_interactive_gui()
    else:
        pass
        # cct.generate_mesh_and_save()
        # cct.load_pro_and_solve()
