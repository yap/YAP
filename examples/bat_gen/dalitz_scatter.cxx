void dalitz_scatter()
{

    TFile* F = TFile::Open("DKSPIPI_mcmc.root", "READ");
    TTree* T = 0;
    F->GetObject("DKSPIPI_mcmc", T);
    if (!T)
        return;

    double m2_01, m2_12;
    unsigned Iteration;
    int Phase;

    T->SetBranchAddress("Phase", &Phase);
    T->SetBranchAddress("m2_01", &m2_01);
    T->SetBranchAddress("m2_12", &m2_12);
    T->SetBranchAddress("Iteration", &Iteration);

    TGraph* g = new TGraph();

    TH2D* h = new TH2D("h", ";m^{2}(#pi^{0}#pi^{0}) [GeV^{2}];m^{2}(K^{0}_{S}#pi^{0}) [GeV^{2}]",
                       100, 0., 1.9, 100, 0.3, 3.05);

    int n = T->GetEntries() / 1265;

    for (long int i = 0; i < T->GetEntries(); ++i) {

        T->GetEntry(i);

        if (Phase <= 0)
            continue;

        double m2_02 = pow(1.86484, 2) + 2 * pow(0.1349766, 2) + pow(0.497614, 2) - m2_01 - m2_12;

        h->Fill(m2_12, m2_01);
        h->Fill(m2_12, m2_02);

        if (Iteration % n != 0)
            continue;

        g->SetPoint(g->GetN(), m2_12, m2_01);

        g->SetPoint(g->GetN(), m2_12, m2_02);

    }

    gStyle->SetPalette(56, 0);
    h->SetStats(false);
    h->Draw("colz");

    g->SetMarkerStyle(7);
    g->SetMarkerColor(4);
    g->Draw("sameP");
    std::cout << g->GetN() << std::endl;
}
