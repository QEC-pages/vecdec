import csv
import pathlib
import time
from dataclasses import dataclass
from typing import Callable, List, Dict, Any, Set, FrozenSet, Iterable, Tuple
import math
import pymatching
import networkx as nx
import stim
import matplotlib.pyplot as plt


# Define some data for working with the three edge orientations.
@dataclass
class EdgeType:
    pauli: str
    hex_to_hex_delta: complex
    hex_to_qubit_delta: complex
EDGE_TYPES = [
    EdgeType(pauli="X", hex_to_hex_delta=2 - 3j, hex_to_qubit_delta=1 - 1j),
    EdgeType(pauli="Y", hex_to_hex_delta=2 + 3j, hex_to_qubit_delta=1 + 1j),
    EdgeType(pauli="Z", hex_to_hex_delta=4, hex_to_qubit_delta=1),
]
EDGES_AROUND_HEX: List[Tuple[complex, complex]] = [
    (-1 - 1j, +1 - 1j),
    (+1 - 1j, +1),
    (+1, +1 + 1j),
    (+1 + 1j, -1 + 1j),
    (-1 + 1j, -1),
    (-1, -1 - 1j),
]


def generate_circuit_cycle(*,
                           q2i: Dict[complex, int],
                           before_parity_measure_2q_depolarization: float,
                           before_round_1q_depolarization: float,
                           before_cycle_1q_depolarization: float,
                           hex_centers: Dict[complex, int],
                           distance: int,
                           detectors: bool) -> stim.Circuit:
    round_circuits = []
    measurement_times: Dict[FrozenSet[int], int] = {}
    current_time = 0
    measurements_per_round: int
    for r in range(3):

        relevant_hexes = [h for h, category in hex_centers.items() if category == r]

        # Find the edges between the relevant hexes, grouped as X/Y/Z.
        edge_groups: Dict[str, List[FrozenSet[complex]]] = {"X": [], "Y": [], "Z": []}
        for h in relevant_hexes:
            for edge_type in EDGE_TYPES:
                q1 = torus(h + edge_type.hex_to_qubit_delta, distance=distance)
                q2 = torus(h + edge_type.hex_to_hex_delta - edge_type.hex_to_qubit_delta, distance=distance)
                edge_groups[edge_type.pauli].append(frozenset([q1, q2]))
        x_qubits = [q2i[q] for pair in edge_groups["X"] for q in sorted_complex(pair)]
        y_qubits = [q2i[q] for pair in edge_groups["Y"] for q in sorted_complex(pair)]

        circuit = stim.Circuit()
        if before_round_1q_depolarization > 0:
            circuit.append_operation("DEPOLARIZE1", sorted(q2i.values()), before_round_1q_depolarization)

        # Make all the parity operations Z basis parities.
        circuit.append_operation("H", x_qubits)
        circuit.append_operation("H_YZ", y_qubits)

        # Turn parity observables into single qubit observables.
        pair_targets = [
            q2i[q]
            for group in edge_groups.values()
            for pair in group
            for q in sorted_complex(pair)
        ]
        if before_parity_measure_2q_depolarization > 0:
            circuit.append_operation("DEPOLARIZE2", pair_targets, before_parity_measure_2q_depolarization)
        circuit.append_operation("CNOT", pair_targets)

        # Measure
        for k in range(0, len(pair_targets), 2):
            edge_key = frozenset([pair_targets[k], pair_targets[k + 1]])
            measurement_times[edge_key] = current_time
            current_time += 1
        circuit.append_operation("M", pair_targets[1::2])

        # Restore qubit bases.
        circuit.append_operation("CNOT", pair_targets)
        circuit.append_operation("H_YZ", y_qubits)
        circuit.append_operation("H", x_qubits)

        # Multiply relevant measurements into the observable.
        included_measurements = []
        for group in edge_groups.values():
            for pair in group:
                a, b = pair
                if a.real == b.real == 1:
                    edge_key = frozenset([q2i[a], q2i[b]])
                    included_measurements.append(stim.target_rec(measurement_times[edge_key] - current_time))
        circuit.append_operation("OBSERVABLE_INCLUDE", included_measurements, 0)

        round_circuits.append(circuit)
    measurements_per_cycle = current_time
    measurements_per_round = measurements_per_cycle // 3

    # Determine which sets of measurements to compare in order to get detection events in the bulk.
    if detectors:
        for r in range(3):
            circuit = stim.Circuit()
            relevant_hexes = [h for h, category in hex_centers.items() if category == (r + 1) % 3]
            end_time = (r + 1) * measurements_per_round
            for h in relevant_hexes:
                record_targets = []
                for a, b in EDGES_AROUND_HEX:
                    q1 = torus(h + a, distance=distance)
                    q2 = torus(h + b, distance=distance)
                    edge_key = frozenset([q2i[q1], q2i[q2]])
                    relative_index = (measurement_times[edge_key] - end_time) % measurements_per_cycle - measurements_per_cycle
                    record_targets.append(stim.target_rec(relative_index))
                    record_targets.append(stim.target_rec(relative_index - measurements_per_cycle))
                circuit.append_operation("DETECTOR", record_targets, [h.real, h.imag, 0])
            circuit.append_operation("SHIFT_COORDS", [], [0, 0, 1])
            round_circuits[r] += circuit

    full_circuit = stim.Circuit()
    if before_cycle_1q_depolarization > 0:
        full_circuit.append_operation("DEPOLARIZE1", sorted(q2i.values()), before_cycle_1q_depolarization)
    full_circuit += round_circuits[0] + round_circuits[1] + round_circuits[2]
    return full_circuit


def generate_circuit(distance: int, cycles: int,
                     before_parity_measure_2q_depolarization: float,
                     before_round_1q_depolarization: float,
                     before_cycle_1q_depolarization: float,
                     start_of_all_noisy_cycles_1q_depolarization: float,
                     ) -> stim.Circuit:

    # Generate and categorize the hexes defining the circuit.
    hex_centers: Dict[complex, int] = {}
    for row in range(3 * distance):
        for col in range(2 * distance):
            center = row * 2j + 2 * col - 1j * (col % 2)
            category = (-row - col % 2) % 3
            hex_centers[torus(center, distance=distance)] = category

    # Find all the qubit positions around the hexes.
    qubit_coordinates: Set[complex] = set()
    for h in hex_centers:
        for edge_type in EDGE_TYPES:
            for sign in [-1, +1]:
                q = h + edge_type.hex_to_qubit_delta * sign
                qubit_coordinates.add(torus(q, distance=distance))

    # Assign integer indices to the qubit positions.
    q2i: Dict[complex, int] = {q: i for i, q in enumerate(sorted_complex(qubit_coordinates))}

    # Generate a circuit performing the parity measurements that are part of each round.
    # Also keep track of the exact order the measurements occur in.
    round_circuit_no_noise_no_detectors = generate_circuit_cycle(
        q2i=q2i,
        before_parity_measure_2q_depolarization=0,
        before_round_1q_depolarization=0,
        before_cycle_1q_depolarization=0,
        hex_centers=hex_centers,
        distance=distance,
        detectors=False,
    )
    round_circuit_no_noise_yes_detectors = generate_circuit_cycle(
        q2i=q2i,
        before_parity_measure_2q_depolarization=0,
        before_round_1q_depolarization=0,
        before_cycle_1q_depolarization=0,
        hex_centers=hex_centers,
        distance=distance,
        detectors=True,
    )
    round_circuit_yes_noise_yes_detectors = generate_circuit_cycle(
        q2i=q2i,
        before_parity_measure_2q_depolarization=before_parity_measure_2q_depolarization,
        before_round_1q_depolarization=before_round_1q_depolarization,
        before_cycle_1q_depolarization=before_cycle_1q_depolarization,
        hex_centers=hex_centers,
        distance=distance,
        detectors=True,
    )

    # Put together the pieces to get a correctable noisy circuit with noiseless time padding
    # (since the time boundaries are not implemented yet).
    full_circuit = stim.Circuit()
    for q, i in q2i.items():
        full_circuit.append_operation("QUBIT_COORDS", [i], [q.real, q.imag])

    # Initialize data qubits along logical observable column into correct basis for observable to be deterministic.
    qubits_along_column = sorted([q for q in qubit_coordinates if q.real == 1], key=lambda v: v.imag)
    initial_bases_along_column = "ZY_ZX_" * distance
    x_initialized = [q2i[q] for q, b in zip(qubits_along_column, initial_bases_along_column) if b == "X"]
    y_initialized = [q2i[q] for q, b in zip(qubits_along_column, initial_bases_along_column) if b == "Y"]
    full_circuit.append_operation("H", x_initialized)
    full_circuit.append_operation("H_YZ", y_initialized)

    full_circuit += (
            round_circuit_no_noise_no_detectors * 2
            + round_circuit_no_noise_yes_detectors * 2
    )
    if start_of_all_noisy_cycles_1q_depolarization > 0:
        full_circuit.append_operation("DEPOLARIZE1",
                                      sorted(q2i.values()),
                                      start_of_all_noisy_cycles_1q_depolarization)
    full_circuit += (
            round_circuit_yes_noise_yes_detectors * cycles
            + round_circuit_no_noise_yes_detectors * 2
            + round_circuit_no_noise_no_detectors * 2
    )

    # Finish circuit with data measurements.
    qubit_coords_to_measure = [q for q, b in zip(qubits_along_column, initial_bases_along_column) if b != "_"]
    qubit_indices_to_measure= [q2i[q] for q in qubit_coords_to_measure]
    order = {q: i for i, q in enumerate(qubit_indices_to_measure)}
    assert cycles % 2 == 0
    full_circuit.append_operation("H_YZ", y_initialized)
    full_circuit.append_operation("H", x_initialized)
    full_circuit.append_operation("M", qubit_indices_to_measure)

    full_circuit.append_operation("OBSERVABLE_INCLUDE",
                                  [stim.target_rec(i - len(qubit_indices_to_measure)) for i in order.values()],
                                  0)

    return full_circuit


def print_2d(values: Dict[complex, Any]):
    assert all(v.real == int(v.real) for v in values)
    assert all(v.imag == int(v.imag) for v in values)
    assert all(v.real >= 0 and v.imag >= 0 for v in values)
    w = int(max((v.real for v in values), default=0) + 1)
    h = int(max((v.imag for v in values), default=0) + 1)
    s = ""
    for y in range(h):
        for x in range(w):
            s += str(values.get(x + y*1j, "_"))
        s += "\n"
    print(s)


def torus(c: complex, *, distance: int) -> complex:
    r = c.real % (distance * 4)
    i = c.imag % (distance * 6)
    return r + i*1j


def sorted_complex(xs: Iterable[complex]) -> List[complex]:
    return sorted(xs, key=lambda v: (v.real, v.imag))


def run_shots_correct_errors_return_num_correct(circuit: stim.Circuit, num_shots: int):
    """Collect statistics on how often logical errors occur when correcting using detections."""
    e = circuit.detector_error_model()
    m = detector_error_model_to_matching(e)

    t0 = time.monotonic()
    detector_samples = circuit.compile_detector_sampler().sample(num_shots, append_observables=True)
    t1 = time.monotonic()

    num_correct = 0
    for sample in detector_samples:
        actual_observable = sample[-1]
        detectors_only = sample.copy()
        detectors_only[-1] = 0
        predicted_observable = m.decode(detectors_only)[0]
        num_correct += actual_observable == predicted_observable
    t2 = time.monotonic()

    # decode_time = t2 - t1
    # sample_time = t1 - t0
    # print("decode", decode_time, "sample", sample_time)

    return num_correct


def detector_error_model_to_matching(model: stim.DetectorErrorModel) -> pymatching.Matching:
    """Convert stim error model into a pymatching graph."""
    det_offset = 0

    def _iter_model(m: stim.DetectorErrorModel, reps: int, callback: Callable[[float, List[int], List[int]], None]):
        nonlocal det_offset
        for _ in range(reps):
            for instruction in m:
                if isinstance(instruction, stim.DemRepeatBlock):
                    _iter_model(instruction.body_copy(), instruction.repeat_count, callback)
                elif isinstance(instruction, stim.DemInstruction):
                    if instruction.type == "error":
                        dets = []
                        frames = []
                        for t in instruction.targets_copy():
                            v = str(t)
                            if v.startswith("D"):
                                dets.append(int(v[1:]) + det_offset)
                            elif v.startswith("L"):
                                frames.append(int(v[1:]))
                            else:
                                raise NotImplementedError()
                        p = instruction.args_copy()[0]
                        callback(p, dets, frames)
                    elif instruction.type == "shift_detectors":
                        det_offset += instruction.targets_copy()[0]
                    elif instruction.type == "detector":
                        pass
                    elif instruction.type == "logical_observable":
                        pass
                    else:
                        raise NotImplementedError()
                else:
                    raise NotImplementedError()

    g = nx.Graph()
    num_detectors = model.num_detectors
    for k in range(num_detectors):
        g.add_node(k)
    g.add_node(num_detectors, is_boundary=True)
    g.add_node(num_detectors + 1)
    for k in range(num_detectors + 1):
        g.add_edge(k, num_detectors + 1, weight=9999999999)

    def handle_error(p: float, dets: List[int], frame_changes: List[int]):
        if p == 0:
            return
        if len(dets) == 1:
            dets.append(num_detectors)
        if len(dets) != 2:
            return  # Just ignore correlated error mechanisms (e.g. Y errors / XX errors)
        g.add_edge(*dets, weight=-math.log(p), qubit_id=frame_changes)

    _iter_model(model, 1, handle_error)

    return pymatching.Matching(g)


def sample_error_rates(*,
                       probabilities: List[float],
                       diameter_factor: List[int],
                       append: bool,
                       path: str,
                       shots: int,
                       noisy_cycles: int,
                       before_parity_measure_2q_depolarization_factor: float,
                       before_round_1q_depolarization_factor: float,
                       before_cycle_1q_depolarization_factor: float,
                       start_of_all_noisy_cycles_1q_depolarization_factor: float):
    if not pathlib.Path(path).exists():
        append = False
    with open(path, "a" if append else "w") as f:
        if not append:
            print("distance,physical_error_rate,num_shots,num_correct", file=f)
        print("diameter_factors", diameter_factor)
        print("probabilities", probabilities)
        print("num_shots", shots)
        for p in probabilities:
            s = f"physical error rate {p}:"
            s = s.rjust(50)
            print(s , end="")
            for d in diameter_factor:
                circuit = generate_circuit(
                    distance=d,
                    cycles=noisy_cycles,
                    before_cycle_1q_depolarization=before_cycle_1q_depolarization_factor*p,
                    before_round_1q_depolarization=before_round_1q_depolarization_factor*p,
                    before_parity_measure_2q_depolarization=before_parity_measure_2q_depolarization_factor*p,
                    start_of_all_noisy_cycles_1q_depolarization=start_of_all_noisy_cycles_1q_depolarization_factor*p,
                )
                num_correct = run_shots_correct_errors_return_num_correct(
                    num_shots=shots,
                    circuit=circuit,
                )
                print(f" {shots - num_correct}", end="")
                print(f"{d},{p},{shots},{num_correct}", file=f, flush=True)
            print()


@dataclass
class DistanceExperimentData:
    num_shots: int = 0
    num_correct: int = 0

    @property
    def logical_error_rate(self) -> float:
        return (self.num_shots - self.num_correct) / self.num_shots


def round_adjustment(error_rate: float, rounds: int) -> float:
    randomize_rate = min(1, 2*error_rate)
    round_randomize_rate = 1 - (1 - randomize_rate)**(1 / rounds)
    round_error_rate = round_randomize_rate / 2
    return round_error_rate




def plot_data(path: str, title: str, rounds_per_shot: int):
    distance_to_noise_to_results: Dict[int, Dict[float, DistanceExperimentData]] = {}
    with open(path, "r") as f:
        for row in csv.DictReader(f):
            distance = int(row["distance"])
            physical_error_rate = float(row["physical_error_rate"])
            d1 = distance_to_noise_to_results.setdefault(distance, {})
            d2 = d1.setdefault(physical_error_rate, DistanceExperimentData())
            d2.num_shots += int(row["num_shots"])
            d2.num_correct += int(row["num_correct"])

    markers = "_ov*sp^<>12348PhH+xXDd|"
    for distance in sorted(distance_to_noise_to_results.keys()):
        group = distance_to_noise_to_results[distance]
        xs = []
        ys = []
        for physical_error_rate in sorted(group.keys()):
            data = group[physical_error_rate]
            xs.append(physical_error_rate)
            ys.append(round_adjustment(data.logical_error_rate, rounds=rounds_per_shot))
        plt.plot(xs, ys, label=f"diameter_scale_factor={distance}", marker=markers[distance])

    plt.legend()
    plt.loglog()

    def f(p):
        if abs(p * 100 - int(p * 100)) < 1e-5:
            return str(int(p * 100)) + "%"
        r = f"{p:.3%}"
        while r and r[-2:] == "0%":
            r = r[:-2] + "%"
        return r
    ticks_y = [k*10**-p for k in range(1, 10) for p in range(1, 5) if k*10**-p <= 0.5]
    ticks_x = [k*10**-p for k in range(1, 10) for p in range(1, 5) if k*10**-p <= 0.5]
    ticks_x.extend([p/100 for p in range(12, 20, 2)])
    plt.xticks([x for x in ticks_x], labels=[f(x) for x in ticks_x], rotation=45)
    plt.yticks([y for y in ticks_y], labels=[f(y) for y in ticks_y])
    plt.ylim(0.0001, 0.5)
    plt.xlim(0.001, 0.5)
    plt.title(title)
    plt.ylabel("Logical Error Rate (Vertical Observable)")
    plt.xlabel("Physical Error Rate Parameter")
    plt.grid()
    plt.show()


def sample_single_depolarizing_layer_circuit():
    sample_error_rates(
        shots=20000,
        probabilities=[
            0.001,
            0.01,
            0.02,
            0.03,
            0.04,
            0.05,
            0.06,
            0.07,
            0.08,
            0.09,
            0.10,
            0.11,
            0.12,
        ],
        before_cycle_1q_depolarization_factor=0,
        before_parity_measure_2q_depolarization_factor=0,
        before_round_1q_depolarization_factor=0,
        noisy_cycles=0,
        start_of_all_noisy_cycles_1q_depolarization_factor=1,
        diameter_factor=[1, 2, 3],
        append=True,
        path="data.csv",
    )


def sample_parity_error_circuit():

    sample_error_rates(
        shots=10000,
        probabilities=[
            0.01,
            0.015,
            0.02,
            0.025,
            0.03,
            0.035,
            0.04,
            0.045,
            0.05,
        ],
        before_cycle_1q_depolarization_factor=0,
        before_parity_measure_2q_depolarization_factor=1,
        before_round_1q_depolarization_factor=0,
        noisy_cycles=6,
        start_of_all_noisy_cycles_1q_depolarization_factor=0,
        diameter_factor=[1, 2],
        append=True,
        path="data_from_parity_errors.csv",
    )


def main():
    # plot_data("data_single.csv",
    #           title="LogLog error rates for toric circuit with single layer of 1q depolarization",
    #           rounds_per_shot=1)
    # return
    # sample_single_depolarizing_layer_circuit()

    sample_parity_error_circuit()

    plot_data("data_from_parity_errors.csv",
              title="LogLog error rates per round for 6 cycle (18 round) toric no-ancilla circuit with 2q depolarization before parity measurements",
              rounds_per_shot=18)


if __name__ == '__main__':
    main()