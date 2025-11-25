from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl
import numpy as np

# Start a Wolfram Engine session
session = WolframLanguageSession()

def heunc(a, q, alpha, gamma, delta, z):
    """Evaluate the HeunC function via Wolfram Engine."""
    return session.evaluate(wl.HeunC(a, q, alpha, gamma, delta, z))

def radial(M,omega,t,r,l):
    m = 2*M
    i = 1j
    if omega**2 - m**2 < 0:
        k = -i*np.sqrt(-(omega**2 - m**2))
    else:
        k = np.sqrt(omega**2 - m**2)
    term7 = -i*M*k
    term8 = i*2*M*(k+np.sqrt(k**2+m**2)) + i*M*(m**2)/k + 1
    term9 = 1
    term10 = 1 + i*4*M*np.sqrt(k**2+m**2)
    term11 = l*(l+1) + 4*M**2*(k**2+(m**2)/2) - i*2*M*(k+np.sqrt(k**2+m**2))
    term12 = 1 - (r/(2*M))
    term0 = np.complex128(heunc(term7, term8, term9, term10, term11, term12))
    term1 = np.exp(-i * np.sqrt(k ** 2 + m ** 2) * t)
    term2 = np.exp(i * np.sqrt(k ** 2 + m ** 2) * r)
    term3 = np.exp(i * k * r)
    term4 = (2 * M - r) ** (4 * i * np.sqrt(k ** 2 + m ** 2) * M)
    term5 = np.exp(2 * i * M * k * (1 - r / M))
    term6 = term1 * term2 * term3 * term4 * term5
    term13 = np.complex128(term6) * term0
    return term13

# Example usage:
print(radial(0.5,10,11,10,0))
# Close session when done
session.terminate()
