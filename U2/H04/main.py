from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QComboBox, QPushButton, QVBoxLayout, QWidget, QMessageBox

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

        # Fórmula de throughput em Mbps
        throughput = n_sc * CPrefix * self.num * n_rb * self.n_cc * mimo * 1e-3 * 0.75 * const
        return throughput
    
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
        mimo_value = self.get_mimo()
        fator = 1 if self.prefix == "normal" else 0.857
        throughput_mcs = throughput_bytable * mimo_value * 1e-3 * self.n_cc * fator

        return throughput_mcs, modulation

class ThroughputApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculadora de Throughput LTE")
        self.init_ui()

    def init_ui(self):
        # Layout principal
        layout = QVBoxLayout()

        # Largura de Banda
        self.bandwidth_label = QLabel("Largura de Banda (MHz):")
        self.bandwidth_input = QComboBox()
        self.bandwidth_input.addItems(["1.4", "3", "5", "10", "15", "20"])
        layout.addWidget(self.bandwidth_label)
        layout.addWidget(self.bandwidth_input)

        # Número de Component Carriers
        self.n_cc_label = QLabel("Número de Component Carriers:")
        self.n_cc_input = QComboBox()
        self.n_cc_input.addItems(["1", "2", "3", "4", "5"])
        layout.addWidget(self.n_cc_label)
        layout.addWidget(self.n_cc_input)

        # Configuração MIMO
        self.mimo_label = QLabel("Configuração MIMO:")
        self.mimo_input = QComboBox()
        self.mimo_input.addItems(["Desabilitado", "2x2", "4x4", "8x8"])
        layout.addWidget(self.mimo_label)
        layout.addWidget(self.mimo_input)

        # Prefixo Cíclico
        self.prefix_label = QLabel("Prefixo Cíclico:")
        self.prefix_input = QComboBox()
        self.prefix_input.addItems(["Normal", "Estendido"])
        layout.addWidget(self.prefix_label)
        layout.addWidget(self.prefix_input)

        # MCS
        self.mcs_label = QLabel("MCS (0-28):")
        self.mcs_input = QComboBox()
        self.mcs_input.addItems(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                                 "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28"])
        layout.addWidget(self.mcs_label)
        layout.addWidget(self.mcs_input)

        # Botão Calcular
        self.calculate_button = QPushButton("Calcular Throughput")
        self.calculate_button.clicked.connect(self.calculate_throughput)
        layout.addWidget(self.calculate_button)

        #Parâmetros de saída

        # Resultado
        self.result_label = QLabel("")
        layout.addWidget(self.result_label)

        # Resultado
        self.result_label_mcs = QLabel("")
        layout.addWidget(self.result_label_mcs)

        # Configuração da janela principal
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def calculate_throughput(self):
        try:
            bandwidth = float(self.bandwidth_input.currentText())
            n_cc = int(self.n_cc_input.currentText())
            mimo = self.mimo_input.currentText()
            prefix = self.prefix_input.currentText().lower()
            mcs = int(self.mcs_input.currentText())

            calculator = LTEThroughputCalculator(bandwidth, n_cc, mimo, prefix, mcs)
            throughput_mcs, modulation = calculator.calculateTput_mcs()
            throughput = calculator.calculate_throughput()

            self.result_label.setText(f"Throughput Teórico: {throughput:.2f} Mbps")
            self.result_label_mcs.setText(f"Throughput TBS: {throughput_mcs:.2f} Mbps | Modulação: {modulation}")
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
