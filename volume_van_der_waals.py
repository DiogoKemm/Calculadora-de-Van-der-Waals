# calculo_vdw.py
import json

CONSTANTES = {
    "N2":    [1.37,  0.0387, 28.014],
    "C4H10": [14.66, 0.1226, 58.124],
}

def calcular_volume(params_json: str) -> str:
    p = json.loads(params_json)

    pressao      = float(p["pressao"])
    grama        = float(p["grama"])
    temp_celsius = float(p["temp_celsius"])
    gas          = p["gas"]
    precisao     = float(p.get("precisao", 1e-6))

    a, b, massa_molar = CONSTANTES[gas]
    mols        = grama / massa_molar
    temp_kelvin = temp_celsius + 273.15
    R           = 0.082057

    def f(V):
        return (pressao + (a * mols**2) / V**2) * (V - mols * b) - mols * R * temp_kelvin

    def dV(V, h=1e-5):
        return (f(V + h) - f(V - h)) / (2 * h)

    # bracket inicial partindo do gás ideal
    gas_ideal = mols * R * temp_kelvin / pressao
    if f(gas_ideal) < 0:
        intervalo_menor, intervalo_maior = gas_ideal, gas_ideal + 1
        while f(intervalo_maior) < 0:
            intervalo_maior += 1
    else:
        intervalo_menor, intervalo_maior = gas_ideal - 1, gas_ideal
        while f(intervalo_menor) > 0:
            intervalo_menor -= 1

    # escolha do ponto inicial (critério de convergência)
    h2      = 1e-3
    d2f_val = (dV(intervalo_menor + h2) - dV(intervalo_menor - h2)) / (2 * h2)
    V0 = intervalo_menor if dV(intervalo_menor) * d2f_val < 0 else intervalo_maior

    # Newton-Raphson
    iteracoes = []
    i, fV = 0, f(V0)
    iteracoes.append({"i": i, "V": V0, "fV": fV})
    while abs(fV) > precisao and i < 30:
        i += 1
        V0 -= fV / dV(V0)
        fV  = f(V0)
        iteracoes.append({"i": i, "V": V0, "fV": fV})

    return json.dumps({
        "volume":       V0,
        "volume_ideal": gas_ideal,
        "mols":         mols,
        "convergiu":    abs(fV) <= precisao,
        "iteracoes":    iteracoes,
    })