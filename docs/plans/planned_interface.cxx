

int main()
{

	FinalStateParticle PiPlus = ParticleFactory::NewFSP(211);
	FinalStateParticle PiMinus = -PiPlus;

	InitialStateParticle* D1 = new InitialStateParticle;
	D1->SetFinalState(PiPlus, PiMinus, PiPlus, PiMinus);

	Resonance * Rho = ParticleFactory::NewResonance(113);
	Rho->SetShape(new BreitWigner);
	Rho->AddChannel(PiPlus, PiMinus);

	D1->AddChannel(Rho, Rho);

	InitialStateParticle* D2 = new InitialStateParticle;
	D2->SetFinalState(PiPlus, PiMinus, PiPlus, KMinus);

	...;

	D2->AddChannel(...);

	// yap::SetInitialState(D1);
	yap::SetInitialState(D2);

	yap::InitialState->FourMomentum(ParticleIndex, DataPoint);
	yap::FourMomentum(ParticleIndex, DataPoint);
	yap::InvariantMass(ParticleCombination, DataPoint);
	{ yap::InitialState->InvariantMass(Particlecombination, DataPoint); }
}
