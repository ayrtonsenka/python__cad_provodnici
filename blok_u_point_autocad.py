import ezdxf
import numpy as np
import pandas as pd
import tkinter as tk
import openpyxl
from pyautocad import Autocad
import time
import tkinter.messagebox
from openpyxl import Workbook
from PIL import Image, ImageTk
import pythoncom
import win32com.client

# IZOHIPSE

def line_segment_intersection_2d(p1, p2, q1, q2, tolerance=1e-6):
    p1 = np.array(p1[:2]); p2 = np.array(p2[:2])
    q1 = np.array(q1[:2]); q2 = np.array(q2[:2])
    u = p2 - p1; v = q2 - q1; w0 = p1 - q1
    a = np.dot(u, u); b = np.dot(u, v); c = np.dot(v, v)
    d = np.dot(u, w0); e = np.dot(v, w0)
    denom = a * c - b * b
    if abs(denom) < tolerance:
        return None
    sc = (b * e - c * d) / denom
    tc = (a * e - b * d) / denom
    if not (0 <= sc <= 1 and 0 <= tc <= 1):
        return None
    point_on_p = p1 + sc * u
    point_on_q = q1 + tc * v
    if np.allclose(point_on_p, point_on_q, atol=tolerance):
        return tuple(point_on_p)
    return None

def compute_station(base_start, base_end, point):
    base_vec = np.array(base_end[:2]) - np.array(base_start[:2])
    point_vec = np.array(point[:2]) - np.array(base_start[:2])
    proj_len = np.dot(point_vec, base_vec) / np.linalg.norm(base_vec)
    return proj_len

def save_points_to_excel(points, filename="preseci.xlsx"):
    wb = Workbook()
    ws = wb.active
    ws.title = "Preseci"
    ws.append(["Station", "Elevation"])
    for pt in points:
        ws.append([pt[0], pt[3]])
    wb.save(filename)
    print(f"\n✅ Station i Elevation sačuvani u Excel fajl: {filename}")

def run_izohipse():
    pythoncom.CoInitialize()
    acad = Autocad(create_if_not_exists=True)
    print("✅ Povezan sa AutoCAD-om.")

    print("👉 Klikni na osnovnu liniju...")
    selection = acad.get_selection("Selektuj osnovnu liniju")

    base_line = None
    for obj in selection:
        if obj.ObjectName == 'AcDbLine':
            base_line = obj
            break

    if not base_line:
        print("⚠ Nije selektovana validna linija.")
        return

    p1 = base_line.StartPoint
    p2 = base_line.EndPoint

    layer_name = input("📝 Unesi ime layera u kom su linije za presek: ").strip()

    doc = acad.app.ActiveDocument
    try:
        doc.SelectionSets.Item("FilteredSet").Delete()
    except:
        pass
    ss = doc.SelectionSets.Add("FilteredSet")
    ss.Select(5)

    relevant_objects = []
    for i in range(ss.Count):
        try:
            obj = ss.Item(i)
            if obj.ObjectName == 'AcDbLine' and obj.Layer == layer_name and obj.Handle != base_line.Handle:
                relevant_objects.append(obj)
        except:
            continue

    points = []
    for obj in relevant_objects:
        q1 = obj.StartPoint; q2 = obj.EndPoint
        intersection_xy = line_segment_intersection_2d(p1, p2, q1, q2)
        if intersection_xy:
            line_vec = np.array(q2[:2]) - np.array(q1[:2])
            point_vec = np.array(intersection_xy) - np.array(q1[:2])
            line_len = np.linalg.norm(line_vec)
            if line_len == 0:
                z_interp = q1[2]
            else:
                t = np.dot(point_vec, line_vec) / (line_len ** 2)
                z_interp = q1[2] + t * (q2[2] - q1[2])
            intersection = (intersection_xy[0], intersection_xy[1], z_interp)
            station = compute_station(p1, p2, intersection)
            points.append((station, intersection[0], intersection[1], intersection[2]))

    if points:
        points.sort(key=lambda pt: pt[0])
        save = input("💾 Sačuvaj u Excel? (da/ne): ").strip().lower()
        if save == "da":
            save_points_to_excel(points)
    else:
        print("⚠ Nema preseka sa objektima u layeru.")

# BLOK U POINT
def layer_exists(doc, layer_name):

    for layer in doc.Layers:
        if layer.Name.lower() == layer_name.lower():
            return True
    return False

def blok_u_point_autocad(layer_name, ime_bloka, attrib_z_name="PNTELEV"):

    pythoncom.CoInitialize()
    acad = win32com.client.Dispatch("AutoCAD.Application")
    doc = acad.ActiveDocument
    msp = doc.ModelSpace

    if not layer_exists(doc, layer_name):
        new_layer = doc.Layers.Add(layer_name)
        new_layer.Color = 1  
    else:
        
        for lyr in doc.Layers:
            if lyr.Name.lower() == layer_name.lower() and lyr.Lock == True:
                lyr.Lock = False

    doc.ActiveLayer = doc.Layers.Item(layer_name)

    vec_vidjeno = []

    for obj in msp:
        if obj.ObjectName == "AcDbBlockReference" and obj.EffectiveName.lower() == ime_bloka.lower():
            x, y, z = obj.InsertionPoint
            z_val = z if z != 0 else None 

            if z_val is None or z_val == 0:
                try:
                    for attrib in obj.GetAttributes():
                        if attrib.TagString.lower() == attrib_z_name.lower():
                            z_val = float(attrib.TextString.replace(',', '.'))
                            break
                except Exception as e:
                    print("⚠ Problem pri čitanju atributa:", e)


            if z_val is None:
                z_val = 0.0

            vec_vidjeno.append((x, y, z_val))

    vec_vidjeno = list(set(vec_vidjeno))

    for koordinata in vec_vidjeno:
        try:
            x, y, z = float(koordinata[0]), float(koordinata[1]), float(koordinata[2])
            safe_array = win32com.client.VARIANT(pythoncom.VT_ARRAY | pythoncom.VT_R8, [x, y, z])
            pt = msp.AddPoint(safe_array)
        except Exception as e:
            print("⚠ Greška pri dodavanju tačke:", (x, y, z), e)

    print(f"✅ Dodato {len(vec_vidjeno)} tačaka u layer '{layer_name}'.")
    doc.Regen(1)


def get_user_inputs():
    result = {"pretvori": None, "offset": None, "layer": None, "ime_bloka": None}

    def on_choice(choice):
        if (offset_entry.get() and layer_entry.get() and block_entry.get() and
            offset_entry.get() != offset_placeholder and
            layer_entry.get() != layer_placeholder and
            block_entry.get() != block_placeholder):
            result["pretvori"] = choice
            result["offset"] = float(offset_entry.get())
            result["layer"] = layer_entry.get()
            result["ime_bloka"] = block_entry.get().strip()
            root.destroy()
        else:
            error_label.config(text="Popunite offset, layer i blok!", fg="red")

    def set_placeholder(entry, placeholder):
        entry.insert(0, placeholder)
        entry.config(fg="gray")
        def on_focus_in(event):
            if entry.get() == placeholder:
                entry.delete(0, "end")
                entry.config(fg="black")
        def on_focus_out(event):
            if entry.get() == "":
                entry.insert(0, placeholder)
                entry.config(fg="gray")
        entry.bind("<FocusIn>", on_focus_in)
        entry.bind("<FocusOut>", on_focus_out)

    root = tk.Tk()
    root.title("Opcije")
    root.geometry("420x350") 

    label = tk.Label(root, text="Da li je potrebno da se blokovi pretvore u tačke?", font=("Arial", 14))
    label.pack(pady=15)

    btn_frame = tk.Frame(root)
    btn_frame.pack(pady=10)

    yes_btn = tk.Button(btn_frame, text="DA", width=15, height=2, bg="#4CAF50", fg="white", font=("Arial", 12, "bold"), command=lambda: on_choice("da"))
    yes_btn.pack(side="left", padx=20)

    no_btn = tk.Button(btn_frame, text="NE", width=15, height=2, bg="#F44336", fg="white", font=("Arial", 12, "bold"), command=lambda: on_choice("ne"))
    no_btn.pack(side="left", padx=20)

    # Offset entry
    offset_label = tk.Label(root, text="Unesi offset:", font=("Arial", 12))
    offset_label.pack(pady=(20, 5))
    offset_entry = tk.Entry(root, font=("Arial", 12), width=20)
    offset_entry.pack()
    offset_placeholder = "npr. 1.0"
    set_placeholder(offset_entry, offset_placeholder)

    # Layer entry
    layer_label = tk.Label(root, text="Unesi ime layer-a:", font=("Arial", 12))
    layer_label.pack(pady=(15, 5))
    layer_entry = tk.Entry(root, font=("Arial", 12), width=20)
    layer_entry.pack()
    layer_placeholder = "npr. points"
    set_placeholder(layer_entry, layer_placeholder)

    # Block entry
    block_label = tk.Label(root, text="Unesi ime bloka:", font=("Arial", 12))
    block_label.pack(pady=(15, 5))
    block_entry = tk.Entry(root, font=("Arial", 12), width=20)
    block_entry.pack()
    block_placeholder = "npr. Pnt"
    set_placeholder(block_entry, block_placeholder)

    error_label = tk.Label(root, text="", font=("Arial", 11))
    error_label.pack(pady=10)

    root.mainloop()
    return result["pretvori"], result["offset"], result["layer"], result["ime_bloka"]


def work(offset=1.0, layer=None):

    def get_polyline_vertices(polyline):
        points = []
        name = polyline.ObjectName
        try:
            coords = list(polyline.Coordinates)
            z = getattr(polyline, "Elevation", 0)
            for i in range(0, len(coords), 2):
                x, y = coords[i], coords[i + 1]
                points.append([x, y, z])
        except Exception as e:
            print(f"⚠ Ne mogu da pročitam koordinate iz {name}: {e}")
        return np.array(points)

    def point_to_polyline_distance(pt, polyline_points):
        min_dist = float("inf"); station = 0; total_len = 0
        for i in range(len(polyline_points) - 1):
            p1 = polyline_points[i]; p2 = polyline_points[i + 1]
            vec = p2 - p1; length = np.linalg.norm(vec[:2])
            if length == 0: continue
            t = np.dot(pt[:2] - p1[:2], vec[:2]) / np.dot(vec[:2], vec[:2])
            t_clamped = max(0, min(1, t))
            closest_point = p1 + t_clamped * vec
            dist = np.linalg.norm(pt[:2] - closest_point[:2])
            if dist < min_dist:
                min_dist = dist
                station = total_len + np.linalg.norm(closest_point[:2] - p1[:2])
            total_len += length
        return min_dist, station

    def collect_points_near_polyline(acad, polyline, offset):
        all_points = []
        polyline_pts = get_polyline_vertices(polyline)
        
        for obj in acad.iter_objects("AcDbPoint"):
            try:
                if obj.Layer.lower() != layer.lower():
                    continue
                coords = list(obj.Coordinates)
                if len(coords) == 3:
                    pt = np.array(coords)
                    dist, station = point_to_polyline_distance(pt, polyline_pts)
                    if dist <= offset:
                        all_points.append((pt, station))
            except Exception as e:
                print("Skipping object due to error:", e)    
        return all_points

    def export_to_excel(points, filename="polyline_points.xlsx"):
        df = pd.DataFrame(
            [(i+1, round(station, 3), round(pt[2], 3)) for i, (pt, station) in enumerate(points)],
            columns=["#Point ID", "Station", "Elevation"]
        )
        df = df.sort_values(by="Station").reset_index(drop=True)
        df.to_excel(filename, index=False, sheet_name="Points")
        print(f"Exported {len(points)} points to {filename}")
        root = tk.Tk(); root.withdraw()
        tk.messagebox.showinfo("Kraj", "Gotovo!")
        root.destroy()

    def main_2():
        pythoncom.CoInitialize()
        try:
            acad = Autocad(create_if_not_exists=True)
            print("AutoCAD Connected:", acad.doc.Name)
        except Exception as e:
            print("Failed to connect to AutoCAD:", e)
            pythoncom.CoUninitialize()
            return

        root = tk.Tk(); root.withdraw()
        tk.messagebox.showinfo("Izbor polilinije", "Izaberite poliliniju u AutoCAD-u i pritisnite Enter.")

        try:
            acad.doc.SelectionSets.Item("PLsel").Delete()
        except: pass

        try:
            sset = acad.doc.SelectionSets.Add("PLsel")
            sset.SelectOnScreen()
        except Exception as e:
            print("Selection error:", e)
            return

        if sset.Count != 1:
            print("Please select exactly one polyline.")
            return

        polyline = sset.Item(0)
        if polyline.ObjectName not in ['AcDbPolyline', 'AcDb2dPolyline', 'AcDb3dPolyline']:
            print(f"Selected object is not a supported polyline (got {polyline.ObjectName})")
            return

        points = collect_points_near_polyline(acad, polyline, offset)
        export_to_excel(points)
        sset.Delete()

    main_2()


def run_blok():
    odgovor, offset, layer, ime_bloka = get_user_inputs()
    if odgovor and odgovor.lower() == 'da':
        blok_u_point_autocad(layer, ime_bloka)  
        time.sleep(10)
        work(offset, layer)
    else:
        work(offset, layer)


def main():
    root = tk.Tk()
    root.title("Izbor režima")
    root.geometry("600x350")

    dvant_label = tk.Label(
        root,
        text="Dvant.rs",
        font=("Arial", 22, "bold"),  
        fg="#33A1E0"  
    )
    dvant_label.place(x=10, y=10)
    
    label = tk.Label(root, text="Da li želite da pokrenete izohipse?", font=("Arial", 14))
    label.pack(pady=20)

    btn_frame = tk.Frame(root)
    btn_frame.pack(pady=20)

    yes_btn = tk.Button(btn_frame, text="DA", width=15, height=2, bg="#4CAF50", fg="white", command=lambda:[root.destroy(), run_izohipse()])
    yes_btn.pack(side="left", padx=20)

    no_btn = tk.Button(btn_frame, text="NE", width=15, height=2, bg="#F44336", fg="white", command=lambda:[root.destroy(), run_blok()])
    no_btn.pack(side="left", padx=20)

    root.mainloop()

if __name__ == "__main__":
    main()
