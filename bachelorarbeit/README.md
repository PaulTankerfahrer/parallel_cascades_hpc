# Bachelorarbeit: agentenbetriebenes HPC-Labor

**Arbeitsfrage (noch nicht final):** Kann ein agentenbetriebenes HPC-Labor ein
bestehendes Minimalmodell eigenständig so weiterentwickeln, dass es
MD-Referenzen näher kommt, ohne seine Einfachheit zu verlieren?

**Stand (30.09.):** in Planung.

Beschlossen ([`../docs/ENTSCHEIDUNGEN.md`](../docs/ENTSCHEIDUNGEN.md), E19–E20):
- Der Agent arbeitet in einem **eigenen Repo** und sieht dieses Repo nicht.
- Er startet vom Modell im Stand von Tag `kursprojekt-code`, mit allen
  bekannten Fehlern. Die Fehler dienen als Prüfliste.
- Er rechnet auf einem minimalen HPC-Cluster (Kubernetes, zu Hause).
- Das Agenten-Repo wird angelegt, wenn der Plan steht.

| Datei | Inhalt |
|---|---|
| [`zusammenfassung_fuer_ba_planung.md`](zusammenfassung_fuer_ba_planung.md) | Kontext für die Planung: Modell, Vorgeschichte, bekannte Fehler als Prüfliste, was der Agent bekommt, was im Paper passiert. **Nicht an den Agenten geben.** |
