from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QComboBox, QPushButton, QVBoxLayout, QWidget, QMessageBox, QHBoxLayout, QGridLayout

class LTEThroughputCalculator:
    def __init__(self, bandwidth, n_cc, mimo, prefix, mcs, table_path = r'C:\Users\pedro\Documents\UFRN\2024.2\Comunicações Móveis\Projetos\HandsOnCM\U2\H04\TBS_data.npz'):
        self.bandwidth = bandwidth  # Largura de banda em MHz
        self.n_cc = n_cc  # Número de Component Carriers
        self.mimo = mimo
        self.mimo_table = {
            "Desabilitado": 1, "2x2": 2, "4x4": 4, "8x8": 8
        }  # Configuração MIMO
        self.prefix = prefix  # Prefixo cíclico (normal/estendido)
        self.mcs = mcs  # Modulation and Coding Scheme
        self.table_path = table_path

        # Tabelas de referência
        self.rb_table = {
            1.4: 6, 3: 15, 5: 25, 10: 50, 15: 75, 20: 100
        }  # Resource Blocks por banda
        self.efficiency_table = {
            0: 0.1523, 1: 0.2344, 2: 0.3770, 3: 0.6016, 4: 0.8770, 5: 1.1758, 6: 1.4766,
            7: 1.9141, 8: 2.4063, 9: 2.7305, 10: 3.3223, 11: 3.9023, 12: 4.5234, 13: 5.1152,
            14: 5.5547, 15: 5.9102
        }  # Eficiência por MCS

    def get_resource_blocks(self):
        """Retorna o número de Resource Blocks para a largura de banda selecionada."""
        return self.rb_table.get(self.bandwidth, 0)

    def get_efficiency(self):
        """Retorna a eficiência espectral com base no MCS."""
        return self.efficiency_table.get(self.mcs, 0)

    def get_mimo(self):
        """Retorna o valor MIMO baseado na configuração."""
        return self.mimo_table.get(self.mimo, 0)

    def calculate_throughput(self):
        """Calcula o throughput teórico baseado nos parâmetros fornecidos."""
        n_rb = self.get_resource_blocks()
        #efficiency = self.get_efficiency()
        mimo = self.get_mimo()
        CPrefix = 7 if self.prefix == "normal" else 6
        symbols = 14 if self.prefix == "normal" else 12
        n_sc = 12  # Subcarriers por Resource Block
        const = 2
        re = CPrefix * n_sc

        # Fórmula de throughput em Mbps
        throughput = n_sc * CPrefix * self.num * n_rb * self.n_cc * mimo * 1e-3 * 0.75 * const
        return throughput, n_rb, re
    
    def calculateTput_mcs(self):
        import numpy as np
        imcs = self.mcs
        if imcs <= 8:
            itbs = imcs
            modulation = "QPSK"
            self.num = 2
        elif imcs == 9:
            itbs == imcs
            modulation = "QPSK"
            self.num = 2
        elif imcs == 10:
            itbs = imcs
            modulation = "16-QAM"
            self.num = 4
        elif imcs > 10 and imcs < 16:
            itbs = imcs - 1
            modulation = "16-QAM"
            self.num = 4
        elif imcs == 16:
            itbs = imcs - 1
            modulation = "16-QAM"
            self.num = 4
        elif imcs == 17:
            itbs = imcs - 2
            modulation = "64-QAM"
            self.num = 6
        elif imcs > 17:
            itbs = imcs - 2
            modulation = "64-QAM"
            self.num = 6

        tbs_data = np.load(self.table_path)
        tbs_matrix = tbs_data['TBS']
        itbs_prb = self.get_resource_blocks()
        throughput_bytable = tbs_matrix[itbs][itbs_prb-1]
        tbs_value = throughput_bytable
        mimo_value = self.get_mimo()
        fator = 1 if self.prefix == "normal" else 0.857
        throughput_mcs = throughput_bytable * mimo_value * 1e-3 * self.n_cc * fator

        return throughput_mcs, modulation, self.num, tbs_value, itbs

class ThroughputApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculadora de Throughput LTE")
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
        self.bandwidth_input.addItems(["1.4", "3", "5", "10", "15", "20"])
        input_grid.addWidget(self.bandwidth_input, 0, 1)

        # Prefixo Cíclico
        input_grid.addWidget(QLabel("CP:"), 0, 2)
        self.cp_input = QComboBox()
        self.cp_input.addItems(["Normal", "Estendido"])
        input_grid.addWidget(self.cp_input, 0, 3)

        # MCS
        input_grid.addWidget(QLabel("MCS:"), 1, 0)
        self.mcs_input = QComboBox()
        self.mcs_input.addItems([str(i) for i in range(29)])
        input_grid.addWidget(self.mcs_input, 1, 1)

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

        # Botão Calcular
        self.calculate_button = QPushButton("CALCULAR")
        self.calculate_button.clicked.connect(self.calculate_throughput)  # Conecta ao cálculo
        input_grid.addWidget(self.calculate_button, 2, 3)

        main_layout.addLayout(input_grid)

        # ====== Parâmetros de Saída ======
        main_layout.addWidget(QLabel("Parâmetros de saída", self))
        output_grid = QGridLayout()

        # Campos de saída
        output_grid.addWidget(QLabel("PRB:"), 0, 0)
        self.prb_output = QLineEdit()
        self.prb_output.setReadOnly(True)
        output_grid.addWidget(self.prb_output, 0, 1)

        output_grid.addWidget(QLabel("TBS INDEX:"), 0, 2)
        self.tbs_index_output = QLineEdit()
        self.tbs_index_output.setReadOnly(True)
        output_grid.addWidget(self.tbs_index_output, 0, 3)

        output_grid.addWidget(QLabel("VALOR TBS:"), 1, 0)
        self.tbs_value_output = QLineEdit()
        self.tbs_value_output.setReadOnly(True)
        output_grid.addWidget(self.tbs_value_output, 1, 1)

        output_grid.addWidget(QLabel("MODULAÇÃO:"), 1, 2)
        self.modulation_output = QLineEdit()
        self.modulation_output.setReadOnly(True)
        output_grid.addWidget(self.modulation_output, 1, 3)

        output_grid.addWidget(QLabel("Nº RE:"), 2, 0)
        self.nr_output = QLineEdit()
        self.nr_output.setReadOnly(True)
        output_grid.addWidget(self.nr_output, 2, 1)

        output_grid.addWidget(QLabel("QTD SÍMBOLOS:"), 2, 2)
        self.symbols_output = QLineEdit()
        self.symbols_output.setReadOnly(True)
        output_grid.addWidget(self.symbols_output, 2, 3)

        main_layout.addLayout(output_grid)

        # Taxa de Transmissão pela equação
        rate_layout_eq = QHBoxLayout()
        rate_layout_eq.addWidget(QLabel("Taxa de Transmissão pela equação:"))
        self.rate_output_eq = QLineEdit()
        self.rate_output_eq.setReadOnly(True)
        rate_layout_eq.addWidget(self.rate_output_eq)
        main_layout.addLayout(rate_layout_eq)

        # Taxa de Transmissão
        rate_layout = QHBoxLayout()
        rate_layout.addWidget(QLabel("Taxa de Transmissão:"))
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
            prefix = self.cp_input.currentText().lower()
            mcs = int(self.mcs_input.currentText())

            calculator = LTEThroughputCalculator(bandwidth, n_cc, mimo, prefix, mcs)
            throughput_mcs, modulation, num, tbs, tbs_index = calculator.calculateTput_mcs()
            throughput, prb, re = calculator.calculate_throughput()

            self.symbols_output.setText(f"{num}")
            self.modulation_output.setText(f"{modulation}")
            self.prb_output.setText(f"{prb}")
            self.tbs_value_output.setText(f"{tbs}")
            self.nr_output.setText(f"{re}")
            self.tbs_index_output.setText(f"{tbs_index}")
            self.rate_output_eq.setText(f"{throughput:.2f} Mbps")
            self.rate_output.setText(f"{throughput_mcs:.2f} Mbps")
            #self.result_label_mcs.setText(f"Throughput TBS: {throughput_mcs:.2f} Mbps | Modulação: {modulation}")
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
