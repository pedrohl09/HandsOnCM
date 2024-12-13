# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 11:37:06 2024

@author: Pedro
"""

from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QComboBox, QPushButton, QVBoxLayout, QWidget, QMessageBox, QHBoxLayout, QGridLayout

class NRThroughputCalculator:
    def __init__(self, bandwidth, n_cc, mimo, fr, modulation, numerology, scs, scaling_factor):
        self.bandwidth = bandwidth  # Largura de banda em MHz
        self.n_cc = n_cc  # Número de Component Carriers
        self.mimo = mimo
        self.scs = scs
        self.scaling_factor = scaling_factor
        self.mimo_table = {
            "Desabilitado": 1, "2x2": 2, "4x4": 4, "8x8": 8
        }  # Configuração MIMO
        self.fr = fr  # Frequency Range (1 or 2)
        self.modulation = modulation  # Modulation and Coding Scheme
        self.numerology = numerology

        # Tabelas de referência
        self.scs_15 = {
            5: 25, 10: 52, 15: 79, 20: 106, 25: 133, 30: 100, 40: 216, 50: 270, 60: "N/A", 80: "N/A", 90: "N/A", 100: "N/A"
        }  # Resource Blocks por banda
        
        self.scs_30 = {
            5: 11, 10: 24, 15: 38, 20: 51, 25: 65, 30: 78, 40: 106, 50: 133, 60: 162, 80: 217, 90: 245, 100: 273
        }
        
        self.scs_60 = {
            5: "N/A", 10: 11, 15: 18, 20: 24, 25: 31, 30: 38, 40: 51, 50: 65, 60: 79, 80: 107, 90: 121, 100: 135
        }

        self.mod_order = {
            "qpsk": 2, "16-qam": 4, "32-qam": 6, "64-qam": 6, "256-qam": 8, "1024-qam": 10
        }

    def get_resource_blocks(self):
        """Retorna o número de Resource Blocks baseado no SCS e largura de banda."""
        scs_table = {15: self.scs_15, 30: self.scs_30, 60: self.scs_60}
        table = scs_table.get(self.scs, {})
        return table.get(self.bandwidth, "N/A")

    def get_mimo(self):
        """Retorna o valor MIMO baseado na configuração."""
        return self.mimo_table.get(self.mimo, 0)
    
    def get_mod_order(self):
        """Retorna o valor MIMO baseado na configuração."""
        return self.mod_order.get(self.modulation, 0)

    def calculate_throughput(self):
        """Calcula o throughput teórico baseado nos parâmetros fornecidos."""
        n_rb = self.get_resource_blocks()
        if n_rb == "N/A":
            raise ValueError("Combinação de SCS e largura de banda inválida.")
        mimo = self.get_mimo()
        oh_dl = 0.14 if self.fr == "fr1" else 0.18
        oh_ul = 0.08 if self.fr == "fr2" else 0.1
        order = self.get_mod_order()
        r = 948/1024
        sf = self.scaling_factor
        num = self.numerology
        Ts = (10e-3)/(14*2**num)
        ca = self.n_cc

        # Fórmula de throughput em Mbps
        throughput_dl = mimo * ca * order * sf * r * ((n_rb*12)/Ts) * (1 - oh_dl) * 10e-6
        throughput_ul = mimo * ca * order * sf * r * ((n_rb*12)/Ts) * (1 - oh_ul) * 10e-6
        return throughput_dl, n_rb, throughput_ul
    

class ThroughputApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculadora de Throughput 5G NR")
        self.setGeometry(100, 100, 600, 400)
        self.init_ui()

    def init_ui(self):
        # Widget central
        container = QWidget()
        self.setCentralWidget(container)

        # Layout principal
        main_layout = QVBoxLayout()
        container.setLayout(main_layout)

        # ====== Parâmetros de Entrada ======
        main_layout.addWidget(QLabel("Parâmetros de entrada", self))
        input_grid = QGridLayout()

        # Largura de Banda
        input_grid.addWidget(QLabel("BW (MHz):"), 0, 0)
        self.bandwidth_input = QComboBox()
        self.bandwidth_input.addItems(["5", "10", "15", "20", "25", "30", "40", "50", "60", "80", "90", "100"])
        input_grid.addWidget(self.bandwidth_input, 0, 1)

        # FR
        input_grid.addWidget(QLabel("Frequency Range:"), 0, 2)
        self.FR_input = QComboBox()
        self.FR_input.addItems(["FR1", "FR2"])
        input_grid.addWidget(self.FR_input, 0, 3)

        # SCS
        input_grid.addWidget(QLabel("SCS (MHz):"), 1, 0)
        self.scs_input = QComboBox()
        self.scs_input.addItems(["15", "30", "60"])
        input_grid.addWidget(self.scs_input, 1, 1)

        # Configuração MIMO
        input_grid.addWidget(QLabel("MIMO:"), 1, 2)
        self.mimo_input = QComboBox()
        self.mimo_input.addItems(["Desabilitado", "2x2", "4x4", "8x8"])
        input_grid.addWidget(self.mimo_input, 1, 3)

        # Número de Component Carriers (CA)
        input_grid.addWidget(QLabel("CA:"), 2, 0)
        self.ca_input = QComboBox()
        self.ca_input.addItems(["1", "2", "3", "4", "5"])
        input_grid.addWidget(self.ca_input, 2, 1)

        # Modulation order
        input_grid.addWidget(QLabel("Modulação:"), 2, 2)
        self.mod_order_input = QComboBox()
        self.mod_order_input.addItems([ "QPSK", "16-QAM", "32-QAM", "64-QAM", "256-QAM", "1024-QAM"])
        input_grid.addWidget(self.mod_order_input, 2, 3)

        # Numerologia
        input_grid.addWidget(QLabel("Numerologia:"), 3, 0)
        self.numerology_input = QComboBox()
        self.numerology_input.addItems(["0", "1", "2", "3", "4", "5"])
        input_grid.addWidget(self.numerology_input, 3, 1)

        # SF
        input_grid.addWidget(QLabel("Scaling factor:"), 3, 2)
        self.sf_input = QComboBox()
        self.sf_input.addItems(["1", "0.8", "0.75", "0.4"])
        input_grid.addWidget(self.sf_input, 3, 3)

        # Botão Calcular
        self.calculate_button = QPushButton("CALCULAR")
        self.calculate_button.clicked.connect(self.calculate_throughput)  # Conecta ao cálculo
        input_grid.addWidget(self.calculate_button, 4, 3)

        main_layout.addLayout(input_grid)

        # ====== Parâmetros de Saída ======
        main_layout.addWidget(QLabel("Parâmetros de saída", self))
        output_grid = QGridLayout()

        # Campos de saída

        main_layout.addLayout(output_grid)

        # Taxa de Transmissão pela equação
        rate_layout_dl = QHBoxLayout()
        rate_layout_dl.addWidget(QLabel("5G NR Throughput (downlink):"))
        self.rate_output_dl = QLineEdit()
        self.rate_output_dl.setReadOnly(True)
        rate_layout_dl.addWidget(self.rate_output_dl)
        main_layout.addLayout(rate_layout_dl)

        # Taxa de Transmissão
        rate_layout = QHBoxLayout()
        rate_layout.addWidget(QLabel("5G NR Throughput (uplink):"))
        self.rate_output = QLineEdit()
        self.rate_output.setReadOnly(True)
        rate_layout.addWidget(self.rate_output)
        main_layout.addLayout(rate_layout)

        # Espaçamento final para alinhamento
        main_layout.addStretch()


    def calculate_throughput(self):
        try:
            bandwidth = float(self.bandwidth_input.currentText())
            n_cc = int(self.ca_input.currentText())
            mimo = self.mimo_input.currentText()
            fr = self.FR_input.currentText().lower()
            modulation = self.mod_order_input.currentText().lower()
            numerology = int(self.numerology_input.currentText())
            scs = int(self.scs_input.currentText())
            scaling_factor = float(self.sf_input.currentText())

            calculator = NRThroughputCalculator(bandwidth, n_cc, mimo, fr, modulation, numerology, scs, scaling_factor)
            throughput_dl = calculator.calculate_throughput()

            throughput_dl, n_rb, throughput_ul = calculator.calculate_throughput()
            self.rate_output_dl.setText(f"{throughput_dl:.2f} Mbps")
            self.rate_output.setText(f"{throughput_ul:.2f} Mbps")

        except ValueError:
            QMessageBox.critical(self, "Erro", "Por favor, insira valores válidos para os campos numéricos.")
        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Erro ao calcular throughput: {e}")

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    window = ThroughputApp()
    window.show()
    sys.exit(app.exec_())
