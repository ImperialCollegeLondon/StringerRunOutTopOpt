import copy
from utilities.python.tools import *
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

setnames = [
    "elset=uels",
    "elset=boundary",
    "elset=disbond",
    "nset=tip",
    "nset=separation",
    "elset=csys",
    "nset=csys",
    "elset=nelx",
    "elset=nely",
    "nset=ls",
    "nset=da"
]

koutstep = [
    "*Matrix generate, stiffness, element by element, elset=disbond\n",
    "*Matrix output, stiffness, format=labels\n",
]


def initialise_ls():
    pass


class NodeData(object):
    def __init__(self, data):
        self.n = len(data)
        self.data = data

    def print(self):
        out = ["*Node\n"]
        out += self.data
        return out


class ElementData(object):
    def __init__(self, type, data, ids=[]):
        self.n = len(data)
        self.type = type
        self.data = data
        self.ids = ids

    def print(self):
        out = ["*Element, type=" + self.type + "\n"]
        out += self.data
        return out


class UserElementData(ElementData):
    def __init__(self, type, nnodes, data, ids, props, setname, useProps=True):
        ElementData.__init__(self, type, data, ids)
        self.props = props
        self.useProps = useProps
        self.nnodes = nnodes
        self.setname = setname

    def printProps(self):
        out = ["*User Element, nodes=" + str(self.nnodes) +
               ", type=" + self.type + ", properties=" + str(len(self.props)) + ", coordinates=3\n"]
        out += ["1, 2, 3, 4, 5, 6, 11, 12\n"]
        out += ["5, 1, 2, 3, 4, 5, 6\n"]
        if self.useProps:
            out += ["*UEL Property, elset=" + self.setname + "\n"]
            out += self.props
        return out


class AbqPart(object):
    def __init__(self, name, node_data, elm_data, data):
        self.name = name
        self.nodes = NodeData(node_data)
        self.elements = []
        for elm_entry in elm_data.keys():
            self.elements += [ElementData(elm_entry, elm_data[elm_entry])]
        self.data = data

    def print(self):
        out = ["*Part, name=" + self.name + "\n"]
        out += self.nodes.print()
        for elm_entry in self.elements:
            out += elm_entry.print()
        out += self.data
        out += ["*End Part" + "\n"]
        return out


class AbqSet(object):
    def __init__(self, name, data):
        self.name = name
        self.data = data


class AbqAssembly(object):
    def __init__(self, name, sets, data):
        self.name = name
        self.sets = sets
        self.data = data

    def print(self):
        out = ["*Assembly, name=" + self.name + "\n"]
        out += self.data
        out += ["*End Assembly" + "\n"]
        return out


class AbqStep(object):
    def __init__(self, name, data):
        self.name = name
        self.data = data

    def print(self):
        out = ["*Step, name=" + self.name + "\n"]
        out += self.data
        out += ["*End Step" + "\n"]
        return out


class AbqInterfaceData(object):
    def __init__(self):
        self.nnodel = [18, 8]
        self.coords = []
        self.ndsg = []
        self.uels = []
        self.uelsg = []
        self.guels = []
        self.guelsg = []
        self.con = []
        self.gcon = []
        self.boundary = []
        self.tip = []
        self.separation = []
        self.props = []
        self.orientations = []
        self.thicknesses = []
        self.transform_shell = []
        self.transform_ls = []
        self.da = 1e-6
        self.nnodesuel = 0
        self.nnodesguel = 0
        self.nnodesls = 0
        self.matels = []
        self.lsnds = []
        self.lsinit = []
        self.adjels = []
        self.adjelsg = []
        self.vizels = []
        self.vizelsg = []
        self.vizelcon = []
        self.vizelnds = []
        self.meshsize = 0
        self.alpha = 0
        self.vfrac = 0

    def print(self):
        lines = []
        lines += ["*nnodel"]
        lines += [str(self.nnodel)[1:-1]]
        lines += ["*coords"]
        lines += [str(len(self.coords))]
        lines += [str(x)[1:-1] for x in self.coords]
        lines += ["*uels"]
        lines += [str(len(self.uels))]
        lines += [str(x) for x in self.uels]
        lines += ["*uelsg"]
        lines += [str(x) for x in self.uelsg]
        if (len(self.guels) > 0):
            lines += ["*guels"]
            lines += [str(len(self.guels))]
            lines += [str(x) for x in self.guels]
            lines += ["*guelsg"]
            lines += [str(x) for x in self.guelsg]
        lines += ["*con"]
        lines += [str(x)[1:-1] for x in self.con]
        if (len(self.guels) > 0):
            lines += ["*gcon"]
            lines += [str(x)[1:-1] for x in self.gcon]
        lines += ["*props"]
        lines += [str(len(self.props))]
        lines += [str(x) for x in self.props]
        lines += ["*orientations"]
        lines += [str(len(self.orientations))]
        lines += [str(x) for x in self.orientations]
        lines += ["*thicknesses"]
        lines += [str(len(self.thicknesses))]
        lines += [str(x) for x in self.thicknesses]
        lines += ["*transform_shell"]
        lines += [str(x)[1:-1] for x in self.transform_shell]
        lines += ["*transform_ls"]
        lines += [str(x)[1:-1] for x in self.transform_ls]
        if (len(self.boundary) > 0):
            lines += ["*boundary"]
            lines += [str(len(self.boundary))]
            lines += [str(x) for x in self.boundary]
        if (len(self.guels) > 0):
            lines += ["*tip"]
            lines += [str(len(self.tip))]
            lines += [str(x) for x in self.tip]
            lines += ["*separation"]
            lines += [str(len(self.separation))]
            lines += [str(x) for x in self.separation]
            lines += ["*da"]
            lines += [str(1e-6)]
        lines += ["*ndsg"]
        lines += [str(x) for x in self.ndsg]
        lines += ["*nnodesuel"]
        lines += [str(self.nnodesuel)]
        lines += ["*nnodesguel"]
        lines += [str(self.nnodesguel)]
        lines += ["*nnodesls"]
        lines += [str(self.nnodesls)]
        if (len(self.guels) > 0):
            lines += ["*matels"]
            lines += [str(len(self.matels))]
            lines += [str(x) for x in self.matels]
        lines += ["*lsnds"]
        lines += [str(len(self.lsnds))]
        lines += [str(x) for x in self.lsnds]
        lines += ["*lsinit"]
        lines += [str(x) for x in self.lsinit]
        if (len(self.guels) > 0):
            lines += ["*adjels"]
            lines += [str(len(self.adjels))]
            lines += [str(x) for x in self.adjels]
            lines += ["*adjelsg"]
            lines += [str(x) for x in self.adjelsg]
        if len(self.vizels) > 0:
            lines += ["*vizels"]
            lines += [str(len(self.vizels))]
            lines += [str(x) for x in self.vizels]
            lines += ["*vizelsg"]
            lines += [str(x) for x in self.vizelsg]
            lines += ["*vizelscon"]
            lines += [str(x)[1:-1] for x in self.vizelcon]
            lines += ["*vizelsnds"]
            lines += [str(len(self.vizelnds))]
            lines += [str(x)[1:-1] for x in self.vizelnds]
        lines += ["*meshsize"]
        lines += [str(self.meshsize)]
        lines += ["*alpha"]
        lines += [str(self.alpha)]
        lines += ["*vfrac"]
        lines += [str(self.vfrac)]
        lines += ["*end"]
        return [x.strip() + "\n" for x in lines]


class AbaqusInp(object):

    # ---
    def __init__(self, filename):
        self.full_name = filename
        self.name = filename.split("_in/")[1]

        with open(filename + ".inp", 'r') as f:
            self.raw = [x for x in f.readlines() if "**" not in x]

        self.header = self.parse_header()
        self.parts, self.partname = self.parse_parts()
        self.assembly, self.sets = self.parse_assembly()
        self.materials, self.props = self.parse_materials()
        self.steps = self.parse_steps()
        self.data = AbqInterfaceData()
        self.separation_map = {}
    # ---

    # ---
    def parse_header(self):
        parts_ind = [i for i, x in enumerate(self.raw) if "*part" in x.lower()]
        return self.raw[:parts_ind[0]]
    # ---

    # ---
    def parse_parts(self):
        parts = {}
        parts_data_aux = [x.strip() for x in self.raw
                          if "*part" in x.lower() or "*end part" in x.lower()]
        parts_ind = [i for i, x in enumerate(self.raw)
                     if "*part" in x.lower() or "*end part" in x.lower()]
        for i_entry in range(int(len(parts_ind)/2)):
            entry_name = parts_data_aux[2*i_entry].split("name=")[1]
            if not "rib" in entry_name:
                raw_aux = self.raw[(parts_ind[2*i_entry]+1):parts_ind[2*i_entry+1]]
                nd_ind = [i for i, x in enumerate(
                    raw_aux) if "*node" in x.lower()]
                el_ind = [i for i, x in enumerate(
                    raw_aux) if "*element" in x.lower()]
                kwd_ind = [i for i, x in enumerate(
                    raw_aux) if "*" in x and "**" not in x and "*node" not in x.lower() and "*element" not in x.lower()]
                nd_data = raw_aux[(nd_ind[0] + 1):el_ind[0]]
                el_data = {}
                for iel_entry in range(len(el_ind)):
                    iel_type = raw_aux[el_ind[iel_entry]].split("type=")[
                        1].strip()
                    if iel_entry + 1 == len(el_ind):
                        iel_data = raw_aux[(el_ind[iel_entry] + 1):kwd_ind[0]]
                    else:
                        iel_data = raw_aux[(
                            el_ind[iel_entry] + 1):el_ind[iel_entry + 1]]
                    el_data[iel_type] = iel_data
                csys_ind = [i for i, x in enumerate(
                    raw_aux) if "*orientation" in x.lower()][0]
                csys_data = [float(x)
                            for x in raw_aux[csys_ind + 1].split(",")]
                orient_ind = [i for i, x in enumerate(
                    raw_aux) if "*shell section, elset=uel_section" in x.lower()]
                orient_data_aux = raw_aux[orient_ind[0]+1:]
                orient_data = []
                thick_data = []
                thick_sum = 0
                for orient_entry in orient_data_aux:
                    line_data_aux = orient_entry.split(",")
                    orient_data += [float(line_data_aux[3])]
                    thick_data += [float(line_data_aux[0])]
                    thick_sum += thick_data[-1]
                thick_data = [x/thick_sum for x in thick_data]
                s_data = raw_aux[kwd_ind[0]:]
                parts[entry_name] = AbqPart(entry_name, nd_data, el_data, s_data)
                parts[entry_name].shellcsys = csys_data
                parts[entry_name].orientations = orient_data
                parts[entry_name].thicknesses = thick_data

        return parts, list(parts.keys())[0]
    # ---

    # ---
    def parse_assembly(self):
        sets = {}
        assemb_ind = [i for i, x in enumerate(
            self.raw) if "*assembly" in x.lower() or "*end assembly" in x.lower() or "*submodel" in x.lower()]
        raw_aux = self.raw[(assemb_ind[0] + 1):assemb_ind[1]]

        kwd_ind = [i for i, x in enumerate(
            raw_aux) if "*elset" in x.lower() or "*nset" in x.lower() or "*end assembly" in x.lower()]

        set_ind = [i for i, x in enumerate(
            [raw_aux[j] for j in kwd_ind]) if any(name in x.lower() for name in setnames)]
        set_data = {}
        for set_entry in setnames[:-1]:
            set_data[set_entry.split("=")[1]] = []

        set_names = []
        for set_entry in set_ind:
            set_kwd = raw_aux[kwd_ind[set_entry]]
            if "elset" in set_kwd.lower():
                set_names += [set_kwd.split("elset=")[1].split(",")[0].strip()]
            elif "nset" in set_kwd.lower():
                set_names += [set_kwd.split("nset=")[1].split(",")[0].strip()]
            if "generate" in set_kwd.lower():
                aux_data = [raw_aux[kwd_ind[set_entry] + 1]]
                aux_data = [int(x) for x in aux_data[0].strip().split(",")]
                set_data[set_names[-1]
                         ] = list(range(aux_data[0], aux_data[1]+1, aux_data[2]))
            else:
                aux_data = raw_aux[(kwd_ind[set_entry] + 1)
                                    :kwd_ind[set_entry + 1]]
                set_data[set_names[-1]] = list_from_str(aux_data, int)

        # r = re.compile(setnames[-1])
        # da_matches = [r.search(x) for x in [raw_aux[j] for j in kwd_ind]]
        # da_ind = [i for i, x in enumerate(da_matches) if x is not None]
        # for da_entry in da_ind:
        #     set_kwd = raw_aux[kwd_ind[da_entry]]
        #     set_names += [set_kwd.split("nset=")[1].split(",")[0].strip()]
        #     if "generate" in set_kwd.lower():
        #         aux_data = raw_aux[(kwd_ind[da_entry] + 1)
        #                             :kwd_ind[da_entry + 1]]
        #         aux_data = [int(x) for x in aux_data.split(",")]
        #         set_data[set_names[-1]
        #                  ] = [list(range(aux_data[0], aux_data[1]+1, aux_data[2]))]
        #     else:
        #         aux_data = raw_aux[(kwd_ind[da_entry] + 1)
        #                             :kwd_ind[da_entry + 1]]
        #         set_data[set_names[-1]] = list_from_str(aux_data, int)

        for item in set_data:
            sets[item] = AbqSet(item, set_data[item])
        assemb_name = self.raw[assemb_ind[0]].split("name=")[1].strip()
        assemb = AbqAssembly(assemb_name, sets, raw_aux)
        return assemb, sets
    # ---

    # ---
    def parse_materials(self):
        assemb_ind = [i for i, x in enumerate(
            self.raw) if "*end assembly" in x.lower()]
        stp_ind = [i for i, x in enumerate(self.raw) if "*step" in x.lower()]
        material_data = self.raw[(
            assemb_ind[0]+1):stp_ind[0]] if stp_ind[0] - 1 > assemb_ind[0] else []
        elast_data = []
        if len(material_data) > 0:
            elast_ind = [
                i for i, x in enumerate(material_data) if "*elastic" in x.lower()]
            elast_data = [float(x)
                          for x in material_data[elast_ind[0] + 1].split(",")]
        return material_data, elast_data
    # ---

    # ---
    def parse_steps(self):
        stp_ind = [i for i, x in enumerate(
            self.raw) if "*step" in x.lower() or "*end step" in x.lower()]
        stp_name = self.raw[stp_ind[0]].split("name=")[1].split(",")[0].strip()
        stp_data = self.raw[(stp_ind[0] + 1):stp_ind[1]]
        return [AbqStep(stp_name, stp_data)]
    # ---

    # ---
    def initialise_levelset(self, nhx, nhy):
        plot = False

        if not (nhx == 0 and nhy == 0):
            elx = self.sets['nelx']
            ely = self.sets['nely']
            nx = len(elx.data)
            ny = len(ely.data)

            ndivx = nx - 2
            ndivy = ny - 2
            npixx = ndivx + 1
            npixy = ndivy + 1

            i = Image.open("utilities/python/1x1s.jpg")
            gi = Image.new('RGB', (i.width*nhx, i.height*nhy))
            for ii in range(nhx):
                for jj in range(nhy):
                    gi.paste(i, (ii*i.width, jj*i.height))

            si = gi.resize((npixx, npixy), Image.BILINEAR)
            l_si = si.convert("RGB")

            x = np.linspace(0, nx, nx+1)
            y = np.linspace(0, ny, ny+1)

            X, Y = np.meshgrid(x, y)
            Z = np.zeros(X.shape) + 1
            vals = np.zeros(np.size(X)) + 1

            valid = 0
            for ii in range(1, npixy + 1):
                for jj in range(1, npixx + 1):
                    (r, g, b) = l_si.getpixel((jj-1, ii-1))
                    val = (r+g+b)/3
                    Z[ii, jj] = 2.0*(val/255.0 - 0.5)
                    valid = int(Y[ii, jj] * (nx + 1) + X[ii, jj])
                    vals[valid] = 2.0*(val/255.0 - 0.5)

            if plot:
                plt.imshow(gi)
                plt.show()

                plt.imshow(si)
                plt.show()

                plt.imshow(l_si)
                plt.show()

                plt.contourf(X, Y, Z, np.linspace(-1, 1, 16))
                plt.contour(X, Y, Z, [0], colors="black")
                plt.gca().set_aspect(1)
                plt.show()

            part = self.parts[self.partname]

            x = list_from_str(part.nodes.data, float, True)
            originid = self.sets['csys_origin'].data[0]
            xptid = self.sets['csys_x'].data[0]
            xyptid = self.sets['csys_xy'].data[0]
            origin = x[originid - 1]
            xpt = x[xptid - 1]
            xypt = x[xyptid - 1]

            v1 = [x - y for (x, y) in zip(xpt, origin)]
            v2 = [x - y for (x, y) in zip(xypt, origin)]
            v3 = cross(v1, v2)
            v2 = cross(v3, v1)
            v1 = normalise(v1)
            v2 = normalise(v2)
            v3 = normalise(v3)

            lsids = self.sets['ls'].data
            for iid, ndid in enumerate(lsids):
                if ndid in self.separation_map:
                    lsids[iid] = self.separation_map[ndid]
            lsnds = [x[id - 1] for id in lsids if not id == originid]
            lsx = [[0.0, 0.0, 0.0]] + lsnds

            for ind, nd in enumerate(lsnds):
                nab = dotp(nd, v3)
                ndplane = [x - nab*y for (x, y) in zip(nd, v3)]
                lsx[ind + 1] = [x - y for (x, y) in zip(ndplane, origin)]

            lscsys = [v1, v2, v3]
            for ind, nd in enumerate(lsx):
                lsx[ind] = matvecmul(lscsys, nd)

            lsxnp = np.array(lsx)
            # lsxnp = lsxnp[lsxnp[:, 1].argsort()]
            sortindx = np.lexsort((lsxnp[:, 2], lsxnp[:, 0], lsxnp[:, 1]))
            lsx = lsxnp[sortindx]
            ls = np.array(lsids)[sortindx]

            if plot:
                xgraph = [x[id-1][0] for id in ls]
                ygraph = [x[id-1][1] for id in ls]
                zgraph = [x[id-1][2] for id in ls]
                fig = plt.figure()
                ax = fig.add_subplot(121, projection="3d")
                ax.scatter(xgraph, ygraph, zgraph, c=list(
                    range(len(xgraph))), marker="o", depthshade=False, cmap=plt.get_cmap("winter"))
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.set_zlabel("z")
                ax = fig.add_subplot(122, projection="3d")
                ax.scatter(xgraph, ygraph, zgraph, c=vals.tolist(),
                           marker="o", depthshade=False, cmap=plt.get_cmap("winter"))
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.set_zlabel("z")
                plt.show()
        else:
            ls = np.array(self.sets['ls'].data)
            vals = np.ones(ls.shape)

        return ls.tolist(), vals.tolist()
    # ---

    # ---
    def get_da_pairs(self):
        coords = copy.deepcopy(list_from_str(self.parts[self.partname].nodes.data, float, True))
        da = self.sets['da'].data
        npairs = len(da) // 2
        da11 = self.sets['da11'].data[0]
        da12 = self.sets['da12'].data[0]
        da21 = self.sets['da21'].data[0]
        da22 = self.sets['da22'].data[0]
        pairs = []
        taken = []
        pairs.append([da11, da21])
        taken.append(da11)
        taken.append(da21)
        pairs.append([da12, da22])
        taken.append(da12)
        taken.append(da22)
        while len(pairs) != npairs:
            dir1 = diffv(coords[da12 - 1], coords[da11 - 1])
            dir2 = diffv(coords[da22 - 1], coords[da21 - 1])

            da1 = sumv(coords[da12 - 1], dir1)
            da2 = sumv(coords[da22 - 1], dir2)

            not_taken = [x for x in da if x not in taken]

            min_dist1 = 1000000000.0
            min_dist2 = 1000000000.0
            min_node1 = not_taken[0]
            min_node2 = not_taken[0]
            for node in not_taken:
                x = coords[node - 1]
                d1 = distv(x, da1)
                if d1 < min_dist1:
                    min_dist1 = d1
                    min_node1 = node

                d2 = distv(x, da2)
                if d2 < min_dist2:
                    min_dist2 = d2
                    min_node2 = node

            da11 = da12
            da21 = da22
            da12 = min_node1
            da22 = min_node2

            pairs.append([da12, da22])
            taken.append(da12)
            taken.append(da22)
        
        return pairs
    # ---

    # ---
    def apply_disbond(self):
        sets = self.sets.keys()
        if ("tip" in sets and
            "disbond" in sets and
            "disbond_upper" in sets and
                "separation" in sets):
            part = self.parts[self.partname]
            coords = list_from_str(part.nodes.data, float, True)
            con = list_from_str(part.elements[0].data, int, True)

            delta = 1e-5
            directions = []
            lengths = []
            coords_da = copy.deepcopy(coords)
            pairs_da = self.get_da_pairs()
            for pair in pairs_da:
                tipnode_id = pair[0]
                facenode_id = pair[1]
                tipnode = coords[tipnode_id - 1]
                frontnode = coords[facenode_id - 1]
                vector = diffv(tipnode, frontnode)
                lengths += [norm(vector)]
                directions += [normalise(vector)]
                coords_da[tipnode_id - 1] = sumv(tipnode, scalev(directions[-1], delta))

            # for ipair in range(1, len(tip)+1):
            #     da = self.sets['da' + str(ipair)].data
            #     tipnode_id = da.index(
            #         next(filter(lambda x: x in tip, da), None))
            #     facenode_id = 1 if tipnode_id == 0 else 0
            #     tipnode = coords[da[tipnode_id] - 1]
            #     frontnode = coords[da[facenode_id] - 1]
            #     vector = [x - y for (x, y) in zip(tipnode, frontnode)]
            #     lengths += [math.sqrt(sum([x*x for x in vector]))]
            #     directions += [[x/lengths[-1] for x in vector]]
            #     coords_da[da[tipnode_id] - 1] = [x + (y*delta)
            #                                      for (x, y) in zip(tipnode, directions[-1])]

            disbond_upper_els = {ix: x for (ix, x) in zip(
                self.sets['disbond_upper'].data, [con[y-1] for y in self.sets['disbond_upper'].data])}
            separation = self.sets['separation'].data
            nnodes = len(coords)
            new_coords = coords
            new_con = con
            separation_ndmap = {}
            for inode in separation:
                nnodes += 1
                new_coords += [coords[inode - 1]]
                coords_da += [coords[inode - 1]]
                nd_con_entry = {ix: disbond_upper_els[ix].index(inode)
                                for ix in disbond_upper_els if inode in disbond_upper_els[ix]}
                for iel in nd_con_entry:
                    new_con[iel - 1][nd_con_entry[iel]] = nnodes
                separation_ndmap[inode] = nnodes
            self.separation_map = separation_ndmap

            og_step = copy.deepcopy(self.steps)
            self.steps = [AbqStep("Step-0", koutstep)]

            neldt = ElementData(
                "SC8R",
                str_from_list(new_con, True),
                range(1, len(new_con) + 1)
            )
            part.elements = [neldt]
            part.nodes.data = str_from_list(coords_da, True)

            # remove submodelling
            assembbkup = copy.deepcopy(self.assembly.data)
            smid = [i for i, x in enumerate(
                self.assembly.data) if "submodel" in x.lower()]
            if smid:
                self.assembly.data.remove(self.assembly.data[smid[0]])
                self.assembly.data.remove(self.assembly.data[smid[0]])

            self.print("_da")

            part.nodes.data = str_from_list(new_coords, True)
            self.print()

            self.steps = og_step
            self.assembly.data = assembbkup
            return True

        log("Skipping disbond step -- missing elsets", 3)
        return False
    # ---

    # ---
    def apply_fnm(self, nholesx, nholesy, alpha, vfrac, nsteps):
        sets = self.sets.keys()
        if ("uels" in sets):
            spart = self.parts[self.partname]

            nd_data = list_from_str(spart.nodes.data, float, True)
            offset = len(nd_data) + 1

            el_data = []
            for el_group in spart.elements:
                el_data += list_from_str(el_group.data, int, True)

            uel_ids = self.sets['uels'].data
            uel_data = [el_data[x-1] for x in uel_ids]
            hgl_data = [x[:8] for x in uel_data]
            nel_ids = [ix for ix in range(
                1, len(el_data)+1) if ix not in uel_ids]
            nel_data = [el_data[x-1] for x in nel_ids]

            edges = []
            uel_edges = []
            iedge = 0
            for uel in uel_data:
                aux = []
                e1 = [uel[1-1], uel[2-1]]
                e2 = [uel[2-1], uel[3-1]]
                e3 = [uel[3-1], uel[4-1]]
                e4 = [uel[4-1], uel[1-1]]
                e5 = [uel[5-1], uel[6-1]]
                e6 = [uel[6-1], uel[7-1]]
                e7 = [uel[7-1], uel[8-1]]
                e8 = [uel[8-1], uel[5-1]]
                e1r = [uel[2-1], uel[1-1]]
                e2r = [uel[3-1], uel[2-1]]
                e3r = [uel[4-1], uel[3-1]]
                e4r = [uel[1-1], uel[4-1]]
                e5r = [uel[6-1], uel[5-1]]
                e6r = [uel[7-1], uel[6-1]]
                e7r = [uel[8-1], uel[7-1]]
                e8r = [uel[5-1], uel[8-1]]
                if e1 not in edges and e1r not in edges:
                    edges += [e1]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e1 in edges:
                        aux += [edges.index(e1)]
                    else:
                        aux += [edges.index(e1r)]

                if e2 not in edges and e2r not in edges:
                    edges += [e2]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e2 in edges:
                        aux += [edges.index(e2)]
                    else:
                        aux += [edges.index(e2r)]

                if e3 not in edges and e3r not in edges:
                    edges += [e3]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e3 in edges:
                        aux += [edges.index(e3)]
                    else:
                        aux += [edges.index(e3r)]

                if e4 not in edges and e4r not in edges:
                    edges += [e4]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e4 in edges:
                        aux += [edges.index(e4)]
                    else:
                        aux += [edges.index(e4r)]

                if e5 not in edges and e5r not in edges:
                    edges += [e5]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e5 in edges:
                        aux += [edges.index(e5)]
                    else:
                        aux += [edges.index(e5r)]

                if e6 not in edges and e6r not in edges:
                    edges += [e6]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e6 in edges:
                        aux += [edges.index(e6)]
                    else:
                        aux += [edges.index(e6r)]

                if e7 not in edges and e7r not in edges:
                    edges += [e7]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e7 in edges:
                        aux += [edges.index(e7)]
                    else:
                        aux += [edges.index(e7r)]

                if e8 not in edges and e8r not in edges:
                    edges += [e8]
                    aux += [iedge]
                    iedge += 1
                else:
                    if e8 in edges:
                        aux += [edges.index(e8)]
                    else:
                        aux += [edges.index(e8r)]

                uel_edges += [aux]

            offset2 = offset + len(edges)
            for iuel, uel in enumerate(uel_edges):
                uel_data[iuel] += [x + offset for x in uel] + \
                    [2*iuel + offset2, 2*iuel + 1 + offset2]

            nfloatnds = len(edges) + 2*len(uel_data)
            for ifltnd in range(nfloatnds):
                nd_data += [[0.0]*3]

            spart.elements = []
            if (len(nel_data) > 0):
                neldt = ElementData(
                    "SC8R",
                    str_from_list(nel_data, True, nel_ids),
                    nel_ids
                )
                spart.elements += [neldt]
            ueldt = UserElementData(
                "U810",
                18,
                str_from_list(uel_data, True, uel_ids),
                uel_ids,
                str_from_list(self.props),
                "uels"
            )
            hg_ids = [x + len(nel_ids) + len(uel_ids)
                      for x in range(1, len(uel_ids) + 1)]
            hgldt = ElementData(
                "C3D8RH",
                str_from_list(hgl_data, True, hg_ids),
                uel_ids
            )
            spart.elements += [ueldt]
            spart.elements += [hgldt]
            spart.nodes.data = str_from_list(nd_data, True)
            self.header += ueldt.printProps()

            guel_ids = self.sets['disbond'].data
            guel_data = [el_data[x-1][:8] for x in guel_ids]
            guel_ids = [x + len(nel_ids) + len(uel_ids) + len(hg_ids)
                        for x in range(1, len(guel_ids) + 1)]

            if (len(guel_data) > 0):
                gueldt = UserElementData(
                    "U888",
                    8,
                    str_from_list(guel_data, True, guel_ids),
                    guel_ids,
                    str_from_list(self.props),
                    "disbond",
                    False
                )
                self.header += gueldt.printProps()
                spart.elements += [gueldt]

            ndsuel = set([x for y in uel_data + guel_data for x in y])
            ndmap = {gx: lx for (gx, lx) in zip(
                ndsuel, range(1, len(ndsuel)+1))}
            elmap = {gx: lx for (gx, lx) in zip(
                uel_ids, range(1, len(uel_ids)+1))}
            gelmap = {gx: lx for (gx, lx) in zip(
                guel_ids, range(1, len(guel_ids)+1))}

            nd_data_local = [nd_data[x-1] for x in ndmap]
            uel_ids_local = [elmap[x] for x in uel_ids]

            uel_data_local = [[ndmap[x] for x in y] for y in uel_data]
            guel_data_local = [[ndmap[x] for x in y] for y in guel_data]

            adjel_ids = [
                x for x in self.sets['uels'].data if x in self.sets['disbond'].data]
            adjel_index = [i+1 for i, x in enumerate(
                self.sets['disbond'].data) if x in adjel_ids]
            adjel_ids = [elmap[x] for x in adjel_ids]

            vizel_ids = [x for x in range(1, len(el_data) + 1)
                         if (not x in self.sets['uels'].data)]
            if len(vizel_ids) > 0:
                vizel_data = [el_data[x - 1] for x in vizel_ids]
                vizel_nds = set([x for y in vizel_data for x in y])
                vizel_ndmap = {gx: lx for (gx, lx) in zip(
                    vizel_nds, range(1, len(vizel_nds) + 1))}
                vizel_nd_data = [nd_data[x - 1] for x in vizel_ndmap]
                vizel_ids_new = [x + len(nel_ids) + len(uel_ids) + len(hg_ids) + len(guel_ids)
                                 for x in range(1, len(vizel_ids) + 1)]
                vizel_map = {gx: lx for (gx, lx) in zip(
                    vizel_ids_new, range(1, len(vizel_ids_new) + 1))}
                vizel_ids_local = [vizel_map[x] for x in vizel_ids_new]
                vizel_data_local = [[vizel_ndmap[x]
                                    for x in y] for y in vizel_data]

                vizeldt = UserElementData(
                    "U889",
                    8,
                    str_from_list(vizel_data, True, vizel_ids_new),
                    vizel_ids_new,
                    str_from_list(self.props),
                    "disbond",
                    False
                )
                self.header += vizeldt.printProps()
                spart.elements += [vizeldt]

            ls = []
            lsval = []
            ls, lsval = self.initialise_levelset(nholesx, nholesy)

            self.data.coords = nd_data_local
            self.data.ndsg = [x for x in ndmap]
            self.data.uels = uel_ids_local
            self.data.uelsg = uel_ids
            self.data.guels = [gelmap[x] for x in guel_ids]
            self.data.guelsg = guel_ids
            self.data.con = uel_data_local
            self.data.gcon = guel_data_local
            self.data.transform_shell, self.data.transform_ls, self.data.meshsize = self.compute_csys()
            self.data.orientations = spart.orientations
            self.data.thicknesses = spart.thicknesses
            print(self.sets.keys())
            print(self.sets['tip'].data)
            if 'boundary' in self.sets.keys():
                self.data.boundary = [elmap[x]
                                      for x in self.sets['boundary'].data]
            if 'tip' in self.sets.keys():
                self.data.tip = [ndmap[x] for x in self.sets['tip'].data]
            if 'separation' in self.sets.keys():
                self.data.separation = [ndmap[x]
                                        for x in self.sets['separation'].data]
            self.data.props = self.props
            self.data.nnodesuel = len(set([x for y in uel_data for x in y]))
            self.data.nnodesguel = len(set([x for y in guel_data for x in y]))
            self.data.nnodesls = len(set([x for y in uel_data for x in y[:4]]))
            if 'disbond' in self.sets.keys():
                self.data.matels = self.sets['disbond'].data
            self.data.lsnds = [ndmap[x] for x in ls]
            self.data.lsinit = lsval
            self.data.adjels = adjel_ids
            self.data.adjelsg = adjel_index
            self.data.alpha = alpha
            self.data.vfrac = vfrac
            if len(vizel_ids) > 0:
                self.data.vizels = vizel_ids_local
                self.data.vizelsg = vizel_ids_new
                self.data.vizelcon = vizel_data_local
                self.data.vizelnds = vizel_nd_data

            og_step = self.steps[0]
            for i in range(2, nsteps - 1 + 2):
                name = "Step-"+str(i)
                self.steps += [AbqStep(name, og_step.data)]

            self.print("_fnm")
            self.print_info()
            return True

        log("Skipping fnm step -- missing element set uels", 3)
        return False
    # ---

    # ---
    def compute_csys(self):
        shellcsys = []
        lscsys = []
        part = self.parts[self.partname]

        v1 = [part.shellcsys[x] for x in [0, 1, 2]]
        v2 = [part.shellcsys[x] for x in [3, 4, 5]]
        v3 = cross(v1, v2)
        v2 = cross(v3, v1)
        v1 = normalise(v1)
        v2 = normalise(v2)
        v3 = normalise(v3)
        shellcsys = transpose([v1, v2, v3])

        origin = self.sets['csys_origin']
        xpt = self.sets['csys_x']
        xypt = self.sets['csys_xy']

        x = list_from_str(part.nodes.data, float, True)
        origin = x[origin.data[0] - 1]
        xpt = x[xpt.data[0] - 1]
        xypt = x[xypt.data[0] - 1]

        # print(origin)
        # print(xpt)
        # print(xypt)

        v1 = [x - y for (x, y) in zip(xpt, origin)]
        meshsize = norm(v1)
        v2 = [x - y for (x, y) in zip(xypt, origin)]

        v3 = cross(v1, v2)
        v2 = cross(v3, v1)
        # print(v1)
        # print(v2)
        # print(v3)
        # print("a")
        v1 = normalise(v1)
        # print("b")
        v2 = normalise(v2)
        # print("c")
        v3 = normalise(v3)
        # print("d")
        lscsys = transpose([v1, v2, v3]) + [origin]

        return shellcsys, lscsys, meshsize
    # ---

    # ---
    def print(self, file_extension=""):
        lines = []
        lines += self.header
        for p in self.parts.values():
            lines += p.print()
        lines += self.assembly.print()
        lines += self.materials
        for s in self.steps:
            lines += s.print()
        with open(self.name + file_extension + ".inp", "w+") as f:
            f.writelines(lines)
    # ---

    # ---
    def print_info(self):
        infofile = self.name + "_fnm.info"
        with open(infofile, "w+") as f:
            f.writelines(self.data.print())
