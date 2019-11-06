from mutual_information import MutualInformation
import os


class TestMU:
    """ Это должны были быть юниттесты по заданию, но они не отрисовывают графики, поэтому это просто тесты """

    def test_slides_mu(self):
        """
        Хитмэп по последовательностям со слайдов
        """
        slide_test = ['ACGAAAGU',
                      'UAGUAAUA',
                      'AGGUGACU',
                      'CGGCAAUG',
                      'GUGGGAAC']
        m = MutualInformation(slide_test)
        m.heatmap('Testing with slide examples')

    def test_trna_heatmap_all(self):
        """
        Хитмэп для всех записей. У меня - эукариоты
        """
        path = os.path.dirname(os.path.abspath(__file__))
        m = MutualInformation(os.path.join(path, "eukaryotic-tRNAs.fa"))
        m.heatmap('Complete heatmap of eukarya')

    def test_trna_subset(self):
        """
        Хитмэп для сэмплов (набора столбцов)
        """
        path = os.path.dirname(os.path.abspath(__file__))
        m = MutualInformation(os.path.join(path, "eukaryotic-tRNAs.fa"))
        m.heatmap('Subset heatmap - first N', columns=5)
        m.heatmap('Subset heatmap - array', columns=[0, 2, 4])

    def test_trna_homo(self):
        """
        Хитмэп для homo sapiens
        """
        path = os.path.dirname(os.path.abspath(__file__))
        m = MutualInformation(os.path.join(path, "eukaryotic-tRNAs.fa"), homo=True)
        m.heatmap('Homo sapiens heatmap')

    def test_trna_related(self):
        """
        Хитмэп для родственных видов
        """
        path = os.path.dirname(os.path.abspath(__file__))
        m = MutualInformation(os.path.join(path, "eukaryotic-tRNAs.fa"), related=True)
        m.heatmap('Related species heatmap')

    def run(self):
        """
        Запуск всех тестов. Если нет ошибок - всё ок
        """
        self.test_slides_mu()
        self.test_trna_homo()
        self.test_trna_related()
        self.test_trna_subset()
        self.test_trna_heatmap_all()
