import sys
import os
from ase.build import add_adsorbate
from ase import Atoms
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton,
    QComboBox, QHBoxLayout, QMessageBox, QTableWidget, QTableWidgetItem,
    QFileDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# --- Savaş Gazları ve Metal Oksitler Listesi ---
WAR_GASES_AND_METAL_OXIDES = {
    "Savaş Gazları": {
        "Sarin (GB)": "CP(C)(O)OP(=O)(O)OC(F)(F)F",
        "Soman (GD)": "CCOP(=O)(O)OC(C)(C)C(C)F",
        "Tabun (GA)": "CN(C)P(=O)(O)OC",
        "Hardal Gazı (Sulfur Mustard)": "ClCCSCCCl",
        "Lewisite": "Cl[As](Cl)Cl",
        "Phosgene": "C(=O)ClCl",
        "Cyanogen Chloride": "N#CCl",
        "Hydrogen Cyanide (HCN)": "C#N",
        "VX Gazı": "CCOP(=S)(OCC)SCC[N+](C)(C)C",
        "Chlorine (Cl2)": "ClCl",
        "Arsine": "[AsH3]",
        "Chloropicrin": "C(Cl)(Cl)(Cl)N=O",
        "BZ (3-Quinuclidinyl benzilate)": "O=C(C1CCCCN1)OC2=CC=CC=C2C",
        "Nitrogen Mustard (HN-1)": "ClCCNCCCl",
        "Nitrogen Mustard (HN-2)": "ClCNCCCl",
        "Cyanogen Bromide": "N#CBr",
    },
    "Metal Oksitler": {
        "ZnO": "O[Zn]",
        "TiO2": "O=[Ti]=O",
        "SnO2": "O=[Sn]=O",
        "Fe2O3": "O=[Fe]O[Fe]=O",
        "CuO": "O[Cu]",
        "WO3": "O=[W](=O)=O",
        "Co3O4": "O=[Co]O[Co]=O",
        "NiO": "O[Ni]",
        "V2O5": "O=[V](=O)O[V](=O)=O"
    }
}

def smiles_to_ase_atoms(smiles: str) -> Atoms:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        positions.append([pos.x, pos.y, pos.z])
    positions = np.array(positions)
    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms

class LargeStructureViewer(QWidget):
    def __init__(self, bg_color="black"):
        super().__init__()
        self.setWindowTitle("Yapı Görüntüleyici")
        self.resize(900, 900)
        self.bg_color = bg_color

        self.atoms = None
        self.highlight_index = None
        self.move_mode = False

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        self.viewer = QWebEngineView()
        self.layout.addWidget(self.viewer)

    def set_structure(self, atoms, highlight_index=None, move_mode=False):
        self.atoms = atoms
        self.highlight_index = highlight_index
        self.move_mode = move_mode
        self.display_structure()

    def display_structure(self):
        if self.atoms is None:
            self.viewer.setHtml("<html><body><h3>Yapı yüklenmedi.</h3></body></html>")
            return

        mol_xyz = self.atoms_to_xyz(self.atoms)
        move_mode_js = "true" if self.move_mode else "false"

        highlight_js = ""
        if self.highlight_index is not None and 0 <= self.highlight_index < len(self.atoms):
            highlight_js = f"""
            viewer.setStyle({{serial: {self.highlight_index+1}}}, {{sphere: {{color: 'red', radius: 0.6}}}});
            """

        # 3Dmol.js HTML & JS (atom drag with disabled camera controls when move_mode true)
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
          <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
          <style>html, body {{ margin: 0; height: 100%; overflow: hidden; background-color: {self.bg_color}; }}</style>
        </head>
        <body>
          <div id="viewer" style="width: 100%; height: 100%;"></div>
          <script>
            let viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "{self.bg_color}"}});  
            let data = `{mol_xyz}`;
            let model = viewer.addModel(data, "xyz");
            viewer.setStyle({{}}, {{stick: {{radius:0.2}}, sphere: {{radius:0.5}}}});
            {highlight_js}

            // Koordinat eksenleri
            viewer.addShape({{type:'arrow', start:{{x:0,y:0,z:0}}, end:{{x:2,y:0,z:0}}, color:'red', radius:0.05}});
            viewer.addLabel('X', {{position:{{x:2,y:0,z:0}}, backgroundColor:'red', fontSize:14}});
            viewer.addShape({{type:'arrow', start:{{x:0,y:0,z:0}}, end:{{x:0,y:2,z:0}}, color:'green', radius:0.05}});
            viewer.addLabel('Y', {{position:{{x:0,y:2,z:0}}, backgroundColor:'green', fontSize:14}});
            viewer.addShape({{type:'arrow', start:{{x:0,y:0,z:0}}, end:{{x:0,y:0,z:2}}, color:'blue', radius:0.05}});
            viewer.addLabel('Z', {{position:{{x:0,y:0,z:2}}, backgroundColor:'blue', fontSize:14}});

            viewer.zoomTo();
            viewer.render();

            let moveAtomMode = {move_mode_js};
            let dragging = false;
            let selectedAtom = null;
            let lastMousePos = null;

            // Kamera kontrol ayarları
            viewer.setClickable(!moveAtomMode);
            viewer.setRotate(!moveAtomMode);
            viewer.setZoom(!moveAtomMode);
            viewer.setPan(!moveAtomMode);

            const canvas = viewer.canvas;

            function getMousePos(evt) {{
                let rect = canvas.getBoundingClientRect();
                return {{
                    x: evt.clientX - rect.left,
                    y: evt.clientY - rect.top
                }};
            }}

            function findClosestAtom(x, y) {{
                let minDist = 15;
                let closest = null;
                let atoms = model.selectedAtoms;
                for(let i=0; i<atoms.length; i++) {{
                    let atom = atoms[i];
                    let screen = viewer.project(atom.pos);
                    let dx = screen.x - x;
                    let dy = screen.y - y;
                    let dist = Math.sqrt(dx*dx + dy*dy);
                    if(dist < minDist) {{
                        minDist = dist;
                        closest = atom;
                    }}
                }}
                return closest;
            }}

            // İmleç stili ayarı
            canvas.style.cursor = moveAtomMode ? "pointer" : "grab";

            canvas.onmousedown = function(evt) {{
                if(!moveAtomMode) return;
                let pos = getMousePos(evt);
                selectedAtom = findClosestAtom(pos.x, pos.y);
                if(selectedAtom) {{
                    dragging = true;
                    lastMousePos = pos;
                    canvas.style.cursor = "grabbing";
                }}
            }};

            canvas.onmouseup = function(evt) {{
                if(!moveAtomMode) return;
                dragging = false;
                selectedAtom = null;
                canvas.style.cursor = "pointer";
            }};

            canvas.onmousemove = function(evt) {{
                if(!moveAtomMode) return;
                if(!dragging || !selectedAtom) return;
                let pos = getMousePos(evt);
                let dx = pos.x - lastMousePos.x;
                let dy = pos.y - lastMousePos.y;
                lastMousePos = pos;

                let scale = 0.01;

                // Yalnızca X ve Y düzleminde sürükle
                selectedAtom.pos.x += dx * scale;
                selectedAtom.pos.y -= dy * scale;

                viewer.render();
            }};
          </script>
        </body>
        </html>
        """

        self.viewer.setHtml(html)

    def atoms_to_xyz(self, atoms):
        lines = [str(len(atoms)), "ASE structure"]
        for a in atoms:
            x, y, z = a.position
            lines.append(f"{a.symbol} {x:.6f} {y:.6f} {z:.6f}")
        return "\n".join(lines)

class AtomViewerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Nero Industry - Quantum Espresso Input Generator - Adsorption Edition")
        self.resize(1200, 900)

        self.current_mol_atoms = None
        self.current_slab_atoms = None
        self.current_atoms = None  # adsorbed system
        self.selected_atom_index = None
        self.freeze_indices = set()
        self.structure_type = None  # "molecule", "slab", "adsorbed_system"

        self.qe_exec_path = None
        self.input_folder = None

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # Üst logo + başlık
        top_layout = QHBoxLayout()
        logo_path = "logo.jpg"
        logo_label = QLabel()
        if os.path.exists(logo_path):
            from PyQt5.QtGui import QPixmap
            pixmap = QPixmap(logo_path).scaled(128,128, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(pixmap)
        else:
            logo_label.setText("Logo yok")
        top_layout.addWidget(logo_label, alignment=Qt.AlignLeft | Qt.AlignTop)

        title_sub_layout = QVBoxLayout()
        self.title_label = QLabel("<b>Advanced R&D - Nero Defence Industries</b>")
        self.title_label.setAlignment(Qt.AlignCenter)
        self.subtitle_label = QLabel("Computational Physics Divison")
        self.subtitle_label.setAlignment(Qt.AlignCenter)

        title_sub_layout.addWidget(self.title_label)
        title_sub_layout.addWidget(self.subtitle_label)

        top_layout.addLayout(title_sub_layout)
        self.layout.addLayout(top_layout)

        # API Key
        self.api_label = QLabel("Materials Project API Key:")
        self.api_input = QLineEdit()
        self.api_input.setEchoMode(QLineEdit.Password)
        self.layout.addWidget(self.api_label)
        self.layout.addWidget(self.api_input)

        # Molekül MP ID
        self.mol_label = QLabel("Molekül MP-ID (örn: mp-1234):")
        self.mol_input = QLineEdit()
        self.load_mol_btn = QPushButton("Molekülü Yükle")
        self.load_mol_btn.clicked.connect(self.load_molecule)
        self.layout.addWidget(self.mol_label)
        self.layout.addWidget(self.mol_input)
        self.layout.addWidget(self.load_mol_btn)

        # Slab MP ID ve tekrar sayısı
        self.slab_label = QLabel("Slab MP-ID (örn: mp-1234):")
        self.slab_input = QLineEdit()
        self.repeat_label = QLabel("Slab tekrar sayısı (örn: 1 1 1):")
        self.repeat_input = QLineEdit("1 1 1")
        self.load_slab_btn = QPushButton("Slabı Yükle")
        self.load_slab_btn.clicked.connect(self.load_slab)
        self.layout.addWidget(self.slab_label)
        self.layout.addWidget(self.slab_input)
        self.layout.addWidget(self.repeat_label)
        self.layout.addWidget(self.repeat_input)
        self.layout.addWidget(self.load_slab_btn)

        # Hazır Molekül Seçimi ComboBox
        self.predefined_mol_label = QLabel("Hazır Molekül Seç (Savaş Gazları & Metal Oksitler):")
        self.predefined_mol_combo = QComboBox()
        self.predefined_mol_combo.addItem("Seçiniz")
        for group in WAR_GASES_AND_METAL_OXIDES:
            self.predefined_mol_combo.addItem(f"--- {group} ---")
            for mol_name in WAR_GASES_AND_METAL_OXIDES[group]:
                self.predefined_mol_combo.addItem(mol_name)
        self.predefined_mol_combo.currentTextChanged.connect(self.predefined_mol_selected)
        self.layout.addWidget(self.predefined_mol_label)
        self.layout.addWidget(self.predefined_mol_combo)

        # Atom tablosu
        self.atom_table = QTableWidget()
        self.atom_table.setColumnCount(6)
        self.atom_table.setHorizontalHeaderLabels(["Index", "Element", "X", "Y", "Z", "Dondur"])
        self.atom_table.cellClicked.connect(self.on_atom_selected)
        self.layout.addWidget(self.atom_table)

        # Seçilen atom bilgisi
        self.selected_atom_label = QLabel("Seçili Atom: Yok")
        self.layout.addWidget(self.selected_atom_label)

        # Atom silme butonu
        self.delete_atom_btn = QPushButton("Seçili Atomu Sil")
        self.delete_atom_btn.setEnabled(False)
        self.delete_atom_btn.clicked.connect(self.delete_selected_atom)
        self.layout.addWidget(self.delete_atom_btn)

        # Freeze toggle
        self.freeze_atom_btn = QPushButton("Seçili Atomu Dondur / Çöz")
        self.freeze_atom_btn.setEnabled(False)
        self.freeze_atom_btn.clicked.connect(self.toggle_freeze_atom)
        self.layout.addWidget(self.freeze_atom_btn)

        # Adsorbe parametreleri
        ads_layout = QHBoxLayout()
        self.height_label = QLabel("Adsorbe yüksekliği (Å):")
        self.height_input = QLineEdit("2.0")
        self.vacuum_label = QLabel("Vakum aralığı (Å):")
        self.vacuum_input = QLineEdit("10.0")
        self.create_ads_btn = QPushButton("Adsorbe Sistemi Oluştur")
        self.create_ads_btn.clicked.connect(self.create_adsorbed_system)
        ads_layout.addWidget(self.height_label)
        ads_layout.addWidget(self.height_input)
        ads_layout.addWidget(self.vacuum_label)
        ads_layout.addWidget(self.vacuum_input)
        ads_layout.addWidget(self.create_ads_btn)
        self.layout.addLayout(ads_layout)

        # Molekül taşıma
        move_mol_layout = QHBoxLayout()
        self.move_x_label = QLabel("Molekül Taşı X (Å):")
        self.move_x_input = QLineEdit("0.0")
        self.move_y_label = QLabel("Molekül Taşı Y (Å):")
        self.move_y_input = QLineEdit("0.0")
        self.move_z_label = QLabel("Molekül Taşı Z (Å):")
        self.move_z_input = QLineEdit("0.0")
        self.move_mol_btn = QPushButton("Molekülü Taşı")
        self.move_mol_btn.clicked.connect(self.move_molecule)
        move_mol_layout.addWidget(self.move_x_label)
        move_mol_layout.addWidget(self.move_x_input)
        move_mol_layout.addWidget(self.move_y_label)
        move_mol_layout.addWidget(self.move_y_input)
        move_mol_layout.addWidget(self.move_z_label)
        move_mol_layout.addWidget(self.move_z_input)
        move_mol_layout.addWidget(self.move_mol_btn)
        self.layout.addLayout(move_mol_layout)

        # Arka plan seçici
        bg_layout = QHBoxLayout()
        self.bg_color_combo = QComboBox()
        self.bg_color_combo.addItems(["black", "white", "gray"])
        self.bg_color_combo.currentIndexChanged.connect(self.update_viewer_background)
        bg_layout.addWidget(QLabel("Arka Plan Rengi:"))
        bg_layout.addWidget(self.bg_color_combo)
        self.layout.addLayout(bg_layout)

        # QE input kaydetme
        self.qe_input_btn = QPushButton("Quantum ESPRESSO Input Dosyalarını Kaydet")
        self.qe_input_btn.clicked.connect(self.generate_qe_inputs)
        self.layout.addWidget(self.qe_input_btn)

        # QE yürütülebilir ve run
        qe_layout = QHBoxLayout()
        self.qe_exec_btn = QPushButton("Quantum ESPRESSO Yürütülebilir Dosya Seç")
        self.qe_exec_btn.clicked.connect(self.select_qe_exec)
        qe_layout.addWidget(self.qe_exec_btn)

        self.qe_run_btn = QPushButton("Quantum ESPRESSO Hesaplamasını Başlat")
        self.qe_run_btn.clicked.connect(self.run_qe_calculation)
        self.qe_run_btn.setEnabled(False)
        qe_layout.addWidget(self.qe_run_btn)
        self.layout.addLayout(qe_layout)

        # Hesaplama çıktı tablosu
        self.qe_output_box = QTableWidget()
        self.qe_output_box.setColumnCount(1)
        self.qe_output_box.setHorizontalHeaderLabels(["Çıktı"])
        self.layout.addWidget(self.qe_output_box, stretch=1)

        # Büyük yapı görüntüleyici pencere
        self.large_viewer = LargeStructureViewer(bg_color=self.bg_color_combo.currentText())
        self.large_viewer.show()

        # Atom taşıma modu butonu
        self.move_mode_btn = QPushButton("Atom Taşıma Modunu Aç/Kapat")
        self.move_mode_btn.setCheckable(True)
        self.move_mode_btn.clicked.connect(self.toggle_move_mode)
        self.layout.addWidget(self.move_mode_btn)

        # Atom taşıma mod durum
        self.move_mode = False

        self.update_atom_table()

    def update_viewer_background(self):
        self.large_viewer.bg_color = self.bg_color_combo.currentText()
        self.refresh_viewer()

    def refresh_viewer(self):
        if self.current_atoms is not None:
            self.large_viewer.set_structure(self.current_atoms, self.selected_atom_index, move_mode=self.move_mode)
        elif self.current_slab_atoms is not None:
            self.large_viewer.set_structure(self.current_slab_atoms, move_mode=self.move_mode)
        elif self.current_mol_atoms is not None:
            self.large_viewer.set_structure(self.current_mol_atoms, move_mode=self.move_mode)
        else:
            self.large_viewer.set_structure(None)

    def load_molecule(self):
        mp_id = self.mol_input.text().strip()
        api_key = self.api_input.text().strip()
        if not mp_id or not api_key:
            QMessageBox.warning(self, "Hata", "Lütfen API anahtarı ve molekül MP-ID girin!")
            return
        try:
            with MPRester(api_key) as m:
                entry = m.get_structure_by_material_id(mp_id)
            ase_atoms = AseAtomsAdaptor.get_atoms(entry)
            self.current_mol_atoms = ase_atoms.copy()
            self.structure_type = "molecule"
            self.current_atoms = None
            self.selected_atom_index = None
            self.freeze_indices.clear()
            self.update_atom_table()
            self.refresh_viewer()
        except Exception as e:
            QMessageBox.critical(self, "Hata", f"Molekül yüklenirken hata: {e}")

    def load_slab(self):
        mp_id = self.slab_input.text().strip()
        api_key = self.api_input.text().strip()
        if not mp_id or not api_key:
            QMessageBox.warning(self, "Hata", "Lütfen API anahtarı ve slab MP-ID ile API anahtarı girin!")
            return
        try:
            with MPRester(api_key) as m:
                slab_struct = m.get_structure_by_material_id(mp_id)
            ase_slab = AseAtomsAdaptor.get_atoms(slab_struct)
            # Tekrar sayısı
            rep_str = self.repeat_input.text().strip()
            reps = tuple(int(x) for x in rep_str.split())
            if len(reps) != 3:
                QMessageBox.warning(self, "Hata", "Tekrar sayısı 3 tamsayı olmalı (örn: 1 1 1)")
                return
            slab_repeated = ase_slab.repeat(reps)
            self.current_slab_atoms = slab_repeated.copy()
            self.structure_type = "slab"
            self.current_atoms = None
            self.selected_atom_index = None
            self.freeze_indices.clear()
            self.update_atom_table()
            self.refresh_viewer()
        except Exception as e:
            QMessageBox.critical(self, "Hata", f"Slab yüklenirken hata: {e}")

    def predefined_mol_selected(self, text):
        if text in WAR_GASES_AND_METAL_OXIDES.get("Savaş Gazları", {}) or text in WAR_GASES_AND_METAL_OXIDES.get("Metal Oksitler", {}):
            for group in WAR_GASES_AND_METAL_OXIDES:
                if text in WAR_GASES_AND_METAL_OXIDES[group]:
                    smiles = WAR_GASES_AND_METAL_OXIDES[group][text]
                    self.current_mol_atoms = smiles_to_ase_atoms(smiles)
                    self.structure_type = "molecule"
                    self.current_atoms = None
                    self.selected_atom_index = None
                    self.freeze_indices.clear()
                    self.update_atom_table()
                    self.refresh_viewer()
                    break

    def update_atom_table(self):
        atoms = self.current_atoms if self.current_atoms is not None else self.current_mol_atoms
        if atoms is None:
            atoms = self.current_slab_atoms

        if atoms is None:
            self.atom_table.setRowCount(0)
            self.selected_atom_label.setText("Seçili Atom: Yok")
            self.delete_atom_btn.setEnabled(False)
            self.freeze_atom_btn.setEnabled(False)
            return

        self.atom_table.setRowCount(len(atoms))
        for i, atom in enumerate(atoms):
            self.atom_table.setItem(i, 0, QTableWidgetItem(str(i)))
            self.atom_table.setItem(i, 1, QTableWidgetItem(atom.symbol))
            self.atom_table.setItem(i, 2, QTableWidgetItem(f"{atom.position[0]:.4f}"))
            self.atom_table.setItem(i, 3, QTableWidgetItem(f"{atom.position[1]:.4f}"))
            self.atom_table.setItem(i, 4, QTableWidgetItem(f"{atom.position[2]:.4f}"))
            frozen_item = QTableWidgetItem("Evet" if i in self.freeze_indices else "Hayır")
            frozen_item.setFlags(Qt.ItemIsEnabled)
            self.atom_table.setItem(i, 5, frozen_item)

        self.atom_table.resizeColumnsToContents()
        self.delete_atom_btn.setEnabled(self.selected_atom_index is not None)
        self.freeze_atom_btn.setEnabled(self.selected_atom_index is not None)

    def on_atom_selected(self, row, column):
        self.selected_atom_index = row
        atoms = self.current_atoms if self.current_atoms is not None else self.current_mol_atoms
        if atoms is not None and 0 <= row < len(atoms):
            atom = atoms[row]
            self.selected_atom_label.setText(f"Seçili Atom: {atom.symbol} (Index: {row})")
        self.delete_atom_btn.setEnabled(True)
        self.freeze_atom_btn.setEnabled(True)
        self.refresh_viewer()

    def delete_selected_atom(self):
        idx = self.selected_atom_index
        if idx is None:
            return
        atoms = self.current_atoms if self.current_atoms is not None else self.current_mol_atoms
        if atoms is None or idx >= len(atoms):
            return
        new_atoms = atoms.copy()
        del new_atoms[idx]
        if self.structure_type == "molecule":
            self.current_mol_atoms = new_atoms
            self.current_atoms = None
        elif self.structure_type == "slab":
            self.current_slab_atoms = new_atoms
            self.current_atoms = None
        else:
            self.current_atoms = new_atoms

        self.selected_atom_index = None
        if idx in self.freeze_indices:
            self.freeze_indices.remove(idx)
        self.update_atom_table()
        self.refresh_viewer()

    def toggle_freeze_atom(self):
        idx = self.selected_atom_index
        if idx is None:
            return
        if idx in self.freeze_indices:
            self.freeze_indices.remove(idx)
        else:
            self.freeze_indices.add(idx)
        self.update_atom_table()

    def create_adsorbed_system(self):
        if self.current_mol_atoms is None or self.current_slab_atoms is None:
            QMessageBox.warning(self, "Uyarı", "Önce molekül ve slab yükleyin!")
            return
        try:
            height = float(self.height_input.text())
            vacuum = float(self.vacuum_input.text())
        except Exception:
            QMessageBox.warning(self, "Hata", "Geçerli bir yüksekliği ve vakumu girin!")
            return

        # Vakumlu hücre oluştur
        slab_copy = self.current_slab_atoms.copy()
        slab_copy.set_pbc([True, True, True])
        cell = slab_copy.get_cell()
        cell[2, 2] = cell[2, 2] + vacuum
        slab_copy.set_cell(cell)

        # Molekülü slab üzerine koy
        mol_copy = self.current_mol_atoms.copy()
        # Molekülü slabın üstüne yerleştir (en yüksek atomun z'si baz alınarak)
        z_max_slab = max(slab_copy.get_positions()[:, 2])
        z_min_mol = min(mol_copy.get_positions()[:, 2])
        delta_z = z_max_slab + height - z_min_mol

        mol_copy.translate([0, 0, delta_z])

        # Slab ve molekülü birleştir
        combined_atoms = slab_copy + mol_copy

        self.current_atoms = combined_atoms
        self.structure_type = "adsorbed_system"
        self.selected_atom_index = None
        self.freeze_indices.clear()
        self.update_atom_table()
        self.refresh_viewer()
        QMessageBox.information(self, "Başarılı", "Adsorbe sistemi oluşturuldu!")

    def move_molecule(self):
        if self.current_mol_atoms is None:
            QMessageBox.warning(self, "Uyarı", "Önce molekül yükleyin!")
            return
        try:
            dx = float(self.move_x_input.text())
            dy = float(self.move_y_input.text())
            dz = float(self.move_z_input.text())
        except Exception:
            QMessageBox.warning(self, "Hata", "Geçerli X, Y, Z değerleri girin!")
            return

        self.current_mol_atoms.translate([dx, dy, dz])
        self.refresh_viewer()
        self.update_atom_table()
        QMessageBox.information(self, "Başarılı", f"Molekül ({dx}, {dy}, {dz}) Å taşındı!")

    def toggle_move_mode(self):
        self.move_mode = not self.move_mode
        self.refresh_viewer()

    def select_qe_exec(self):
        exe_path, _ = QFileDialog.getOpenFileName(self, "Quantum ESPRESSO Yürütülebilir Dosya Seç", "", "Executable Files (*)")
        if exe_path:
            self.qe_exec_path = exe_path
            self.qe_run_btn.setEnabled(True)
            QMessageBox.information(self, "Seçildi", f"QE yürütülebilir dosya seçildi:\n{exe_path}")

    def run_qe_calculation(self):
        import subprocess

        if self.qe_exec_path is None:
            QMessageBox.warning(self, "Hata", "Önce QE yürütülebilir dosya seçin!")
            return

        folder = QFileDialog.getExistingDirectory(self, "Input Dosyalarının Bulunduğu Klasörü Seçin")
        if not folder:
            return

        # En azından adsorbe sistemin inputu olmalı
        ads_input = os.path.join(folder, "qe_adsorbed.in")
        if not os.path.exists(ads_input):
            QMessageBox.warning(self, "Hata", "qe_adsorbed.in dosyası seçilen klasörde bulunamadı!")
            return

        # pw.x komutu
        cmd = [self.qe_exec_path, "-in", ads_input]

        try:
            self.qe_output_box.clearContents()
            self.qe_output_box.setRowCount(0)

            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = proc.communicate()

            # Çıktı satırlarını tabloya yaz
            lines = stdout.splitlines()
            for i, line in enumerate(lines):
                self.qe_output_box.insertRow(i)
                self.qe_output_box.setItem(i, 0, QTableWidgetItem(line))

            if proc.returncode != 0:
                QMessageBox.warning(self, "Hata", f"Quantum ESPRESSO çalıştırılırken hata:\n{stderr}")
            else:
                QMessageBox.information(self, "Başarılı", "Quantum ESPRESSO hesaplaması tamamlandı!")

        except Exception as e:
            QMessageBox.critical(self, "Hata", f"Hesaplama sırasında hata: {e}")

    def generate_qe_inputs(self):
        if self.current_atoms is None and self.current_mol_atoms is None and self.current_slab_atoms is None:
            QMessageBox.warning(self, "Uyarı", "Yapı yüklenmedi!")
            return
        folder = QFileDialog.getExistingDirectory(self, "Input Dosyalarının Kaydedileceği Klasörü Seçin")
        if not folder:
            return

        self.input_folder = folder

        def write_input(atoms, filename):
            try:
                cell = atoms.get_cell()
                positions = atoms.get_positions()
                symbols = atoms.get_chemical_symbols()

                # Hücre parametreleri (3 satır)
                cell_params = ""
                for v in cell:
                    cell_params += f"  {v[0]:.8f}  {v[1]:.8f}  {v[2]:.8f}\n"

                # ATOMIC_POSITIONS kısmı
                atomic_positions = ""
                for s, pos in zip(symbols, positions):
                    atomic_positions += f"  {s}  {pos[0]:.8f}  {pos[1]:.8f}  {pos[2]:.8f}\n"

                nat = len(symbols)
                ntyp = len(set(symbols))

                input_text = f"""&control
  calculation = 'scf',
  prefix = '{filename[:-3]}',
  outdir = './tmp/',
  pseudo_dir = './pseudo/',
/
&system
  ibrav = 0,
  nat = {nat},
  ntyp = {ntyp},
  ecutwfc = 40,
  ecutrho = 320,
/
&electrons
  conv_thr = 1.0d-6,
  mixing_beta = 0.7,
/
CELL_PARAMETERS angstrom
{cell_params}
ATOMIC_POSITIONS angstrom
{atomic_positions}
K_POINTS automatic
4 4 4 1 1 1
"""
                with open(os.path.join(folder, filename), "w") as f:
                    f.write(input_text)
            except Exception as e:
                QMessageBox.critical(self, "Hata", f"{filename} dosyası oluşturulurken hata: {e}")

        # Molekül input dosyası
        if self.current_mol_atoms is not None:
            write_input(self.current_mol_atoms, "qe_molecule.in")

        # Slab input dosyası
        if self.current_slab_atoms is not None:
            write_input(self.current_slab_atoms, "qe_slab.in")

        # Adsorbe sistem input dosyası (combined current_atoms)
        if self.current_atoms is not None and self.structure_type == "adsorbed_system":
            write_input(self.current_atoms, "qe_adsorbed.in")

        QMessageBox.information(self, "Başarılı", f"Quantum ESPRESSO input dosyaları oluşturuldu ve klasöre kaydedildi:\n{folder}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AtomViewerApp()
    window.show()
    sys.exit(app.exec_())
