import os
import trimesh
import subprocess

base_path_hd = "/Users/lihongbo/Downloads/research/tmp/"
lp_exe = "/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/src/src"
name_ext = "flipy_emb"

# for obj_id in obj_list:
obj_id = 80557
result_path = "/Users/lihongbo/Downloads/research/tmp/80557/result"
param = ["", "", "", "", "", ""]
input_path = os.path.join(base_path_hd,
                          str(obj_id),
                          "{}_{}.obj".format(obj_id, name_ext))
obj_path = os.path.join(base_path_hd,
                        str(obj_id),
                        "{}_{}.obj".format(obj_id, "sf_norm"))
output_path = os.path.join(result_path,
                           "{}_{}_tri.obj".format(obj_id, name_ext))

iter_num = 400
mesh = trimesh.load(obj_path)
verts_num = len(mesh.vertices)
param[0] = input_path
param[1] = output_path
param[2] = str(8)
param[3] = str(iter_num)
param[4] = str(verts_num)
param[5] = "false"
print(param)
command = [lp_exe]+param
try:
    returncode = subprocess.call(command)
    print("Return code:", returncode)
except Exception as e:
    print("An error occurred:", str(e))
